"""Pure-Python core for the Latch RNA-seq DE workflow.

Everything here is plain pandas / PyDESeq2 so it runs and is tested locally;
the Latch layer in ``wf/__init__.py`` only moves files in and out.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class QCReport:
    n_samples: int
    n_genes_in: int
    n_genes_kept: int
    library_sizes: dict[str, int]
    flagged_samples: list[str]


def load_inputs(counts_path: Path, samples_path: Path, condition_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Counts: genes x samples (first column = gene id). Samples: one row per sample."""
    counts = pd.read_csv(counts_path, index_col=0)
    samples = pd.read_csv(samples_path, index_col=0)
    missing = set(samples.index) - set(counts.columns)
    if missing:
        raise ValueError(f"samples missing from counts matrix: {sorted(missing)}")
    if condition_col not in samples.columns:
        raise ValueError(f"condition column '{condition_col}' not in sample sheet")
    if samples[condition_col].nunique() != 2:
        raise ValueError("this workflow compares exactly two conditions")
    counts = counts[samples.index]
    if (counts < 0).any().any() or not np.allclose(counts, counts.round()):
        raise ValueError("counts must be non-negative integers (raw featureCounts output)")
    return counts.astype(int), samples


def qc_filter(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 2,
              lib_size_z: float = 3.0) -> tuple[pd.DataFrame, QCReport]:
    """Drop low-count genes; flag samples whose library size is an outlier (robust z)."""
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    filtered = counts.loc[keep]
    lib = counts.sum(axis=0)
    log_lib = np.log10(lib.clip(lower=1))
    mad = np.median(np.abs(log_lib - np.median(log_lib))) or 1e-9
    robust_z = 0.6745 * (log_lib - np.median(log_lib)) / mad
    flagged = list(lib.index[np.abs(robust_z) > lib_size_z])
    report = QCReport(
        n_samples=counts.shape[1],
        n_genes_in=counts.shape[0],
        n_genes_kept=int(keep.sum()),
        library_sizes={k: int(v) for k, v in lib.items()},
        flagged_samples=flagged,
    )
    return filtered, report


def run_deseq2(counts: pd.DataFrame, samples: pd.DataFrame, condition_col: str,
               reference: str, alpha: float = 0.05, n_cpus: int = 1) -> pd.DataFrame:
    """PyDESeq2 Wald test, contrast = other level vs ``reference``."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    levels = list(samples[condition_col].unique())
    if reference not in levels:
        raise ValueError(f"reference '{reference}' not in {levels}")
    treated = next(lv for lv in levels if lv != reference)

    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        counts=counts.T,  # PyDESeq2 wants samples x genes
        metadata=samples,
        design=f"~{condition_col}",
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=[condition_col, treated, reference],
                       alpha=alpha, inference=inference, quiet=True)
    stats.summary()
    res = stats.results_df.copy()
    res["significant"] = (res["padj"] < alpha).fillna(False)
    return res.sort_values("padj")


def summarize(res: pd.DataFrame, lfc_threshold: float = 1.0) -> dict:
    sig = res[res["significant"]]
    return {
        "n_tested": int(res["padj"].notna().sum()),
        "n_significant": int(len(sig)),
        "n_up": int((sig["log2FoldChange"] >= lfc_threshold).sum()),
        "n_down": int((sig["log2FoldChange"] <= -lfc_threshold).sum()),
        "top_genes": sig.head(10).index.tolist(),
    }


def volcano_plot(res: pd.DataFrame, out_png: Path, title: str = "") -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = res.dropna(subset=["padj"])
    y = -np.log10(df["padj"].clip(lower=1e-300))
    colors = np.where(df["significant"], "#c0392b", "#95a5a6")
    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)
    ax.scatter(df["log2FoldChange"], y, s=4, c=colors, linewidths=0)
    ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c="k")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title(title or "Differential expression")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    return out_png
