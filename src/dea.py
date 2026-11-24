"""Differential expression utilities for fold-specific feature selection."""

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def run_differential_expression(
    counts: pd.DataFrame,
    conditions: pd.Series,
    alpha: float = 0.05,
    n_cpus: int = 4,
) -> tuple[pd.DataFrame, pd.Index]:
    """Fit pyDESeq2 and return its results and adjusted-p-value-selected genes.

    Counts must contain raw integer counts with samples in rows and genes in
    columns. Conditions must provide matching sample labels for Cancer and
    Healthy. They are aligned to the count rows before fitting.
    """
    if not counts.index.is_unique or not conditions.index.is_unique:
        raise ValueError("Sample identifiers must be unique in counts and conditions.")
    if not counts.columns.is_unique:
        raise ValueError("Gene identifiers must be unique in counts.")
    if set(counts.index) != set(conditions.index):
        raise ValueError("Counts and conditions must contain the same samples.")
    if conditions.isna().any():
        raise ValueError("Conditions must not contain missing values.")
    if set(conditions.unique()) != {"Healthy", "Cancer"}:
        raise ValueError("Conditions must contain both Healthy and Cancer labels only.")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1.")
    if n_cpus < 1:
        raise ValueError("n_cpus must be at least 1.")

    metadata = pd.DataFrame(
        {"Condition": conditions.loc[counts.index].to_numpy()},
        index=counts.index,
    )
    metadata["Condition"] = pd.Categorical(
        metadata["Condition"],
        categories=["Healthy", "Cancer"],
    )

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~Condition",
        refit_cooks=True,
        n_cpus=n_cpus,
    )
    dds.deseq2()

    stats = DeseqStats(
        dds,
        contrast=["Condition", "Cancer", "Healthy"],
        lfc_null=1.0,
        alt_hypothesis="greaterAbs",
    )
    stats.summary()
    stats.lfc_shrink(coeff="Condition[T.Cancer]")
    results = stats.results_df.sort_values("padj", na_position="last").copy()
    significant = results["padj"].notna() & (results["padj"] < alpha)
    selected_genes = results.index[significant]

    return results, selected_genes
