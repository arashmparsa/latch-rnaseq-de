"""Latch workflow: bulk RNA-seq differential expression (PyDESeq2).

Takes the gene x sample counts matrix produced by featureCounts (e.g. from
github.com/arashmparsa/bulk-rnaseq-differential-expression) plus a sample sheet,
and returns a DE table, a QC summary and a volcano plot in Latch Data.
"""
import json
from pathlib import Path

from latch.resources.tasks import medium_task, small_task
from latch.resources.workflow import workflow
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter

from wf import core

metadata = LatchMetadata(
    display_name="Bulk RNA-seq Differential Expression (PyDESeq2)",
    author=LatchAuthor(name="Arash Mehrabi Parsa", github="https://github.com/arashmparsa"),
    license="MIT",
    parameters={
        "counts": LatchParameter(display_name="Counts matrix (CSV, genes x samples)", batch_table_column=True),
        "samples": LatchParameter(display_name="Sample sheet (CSV, sample id first column)"),
        "condition_col": LatchParameter(display_name="Condition column"),
        "reference": LatchParameter(display_name="Reference level (e.g. normal)"),
        "alpha": LatchParameter(display_name="FDR threshold"),
        "output_dir": LatchParameter(display_name="Output directory"),
    },
)


@small_task
def qc_task(counts: LatchFile, samples: LatchFile, condition_col: str) -> LatchFile:
    """Validate inputs, filter low-count genes, flag outlier libraries."""
    df, sheet = core.load_inputs(Path(counts.local_path), Path(samples.local_path), condition_col)
    filtered, report = core.qc_filter(df)
    out = Path("/root/filtered_counts.csv")
    filtered.to_csv(out)
    Path("/root/qc.json").write_text(json.dumps(report.__dict__, indent=2))
    return LatchFile(str(out), "latch:///rnaseq_de_tmp/filtered_counts.csv")


@medium_task
def deseq2_task(filtered_counts: LatchFile, samples: LatchFile, condition_col: str,
                reference: str, alpha: float, output_dir: LatchOutputDir) -> LatchOutputDir:
    """Run PyDESeq2, write results table, summary and volcano plot."""
    import pandas as pd

    counts = pd.read_csv(filtered_counts.local_path, index_col=0)
    sheet = pd.read_csv(samples.local_path, index_col=0)
    res = core.run_deseq2(counts, sheet, condition_col, reference, alpha=alpha, n_cpus=4)

    out = Path("/root/outputs")
    out.mkdir(exist_ok=True)
    res.to_csv(out / "deseq2_results.csv")
    (out / "summary.json").write_text(json.dumps(core.summarize(res), indent=2))
    core.volcano_plot(res, out / "volcano.png", title=f"{condition_col}: vs {reference}")
    return LatchOutputDir(str(out), output_dir.remote_path)


@workflow(metadata)
def rnaseq_de(counts: LatchFile, samples: LatchFile, output_dir: LatchOutputDir,
              condition_col: str = "condition", reference: str = "normal",
              alpha: float = 0.05) -> LatchOutputDir:
    """Bulk RNA-seq differential expression

    # Bulk RNA-seq DE
    Filters low-count genes, flags outlier libraries, then runs a PyDESeq2
    Wald test (treated vs reference) and writes results, summary and a volcano plot.
    """
    filtered = qc_task(counts=counts, samples=samples, condition_col=condition_col)
    return deseq2_task(filtered_counts=filtered, samples=samples, condition_col=condition_col,
                       reference=reference, alpha=alpha, output_dir=output_dir)
