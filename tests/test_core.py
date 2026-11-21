"""Tests on simulated counts with genes we planted as differentially expressed,
so we can measure whether the workflow recovers them."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from wf import core


def simulate(n_genes=2000, n_per_group=4, n_de=100, fold=4.0, seed=0):
    rng = np.random.default_rng(seed)
    base = rng.lognormal(mean=4, sigma=1.2, size=n_genes)
    disp = 0.1
    samples = [f"N{i}" for i in range(n_per_group)] + [f"T{i}" for i in range(n_per_group)]
    cond = ["normal"] * n_per_group + ["tumor"] * n_per_group
    mu = np.tile(base[:, None], (1, len(samples)))
    de_idx = rng.choice(n_genes, n_de, replace=False)
    direction = rng.choice([1, -1], n_de)
    for j, c in enumerate(cond):
        if c == "tumor":
            mu[de_idx, j] *= fold ** direction
    # negative binomial via gamma-Poisson
    lam = rng.gamma(shape=1 / disp, scale=mu * disp)
    counts = rng.poisson(lam)
    genes = [f"g{i}" for i in range(n_genes)]
    counts_df = pd.DataFrame(counts, index=genes, columns=samples)
    sheet = pd.DataFrame({"condition": cond}, index=samples)
    truth = set(np.array(genes)[de_idx])
    return counts_df, sheet, truth


@pytest.fixture(scope="module")
def sim():
    return simulate()


def test_load_inputs_validates(tmp_path):
    counts, sheet, _ = simulate(n_genes=50, n_de=5)
    counts.to_csv(tmp_path / "c.csv")
    sheet.to_csv(tmp_path / "s.csv")
    c, s = core.load_inputs(tmp_path / "c.csv", tmp_path / "s.csv", "condition")
    assert c.shape == (50, 8)
    with pytest.raises(ValueError):
        core.load_inputs(tmp_path / "c.csv", tmp_path / "s.csv", "batch")


def test_qc_flags_broken_library(sim):
    counts, _, _ = sim
    broken = counts.copy()
    broken["T3"] = (broken["T3"] * 0.01).round().astype(int)  # a library that failed
    _, report = core.qc_filter(broken)
    assert report.flagged_samples == ["T3"]


def test_deseq2_recovers_planted_genes(sim):
    counts, sheet, truth = sim
    filtered, _ = core.qc_filter(counts)
    res = core.run_deseq2(filtered, sheet, "condition", reference="normal")
    called = set(res.index[res["significant"]])
    tested_truth = truth & set(filtered.index)
    recall = len(called & tested_truth) / len(tested_truth)
    precision = len(called & tested_truth) / max(len(called), 1)
    print(f"recall={recall:.2f} precision={precision:.2f} called={len(called)}")
    assert recall > 0.8
    assert precision > 0.9


def test_volcano_writes_png(sim, tmp_path):
    counts, sheet, _ = sim
    filtered, _ = core.qc_filter(counts.iloc[:300])
    res = core.run_deseq2(filtered, sheet, "condition", reference="normal")
    out = core.volcano_plot(res, tmp_path / "v.png")
    assert out.exists() and out.stat().st_size > 1000
