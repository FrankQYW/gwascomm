"""Standalone DataFrame-returning CommGWAS enrichment functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ._core import (
    CELL_PHONE_COLUMNS,
    RECEPTOR_COLUMNS,
    LIGAND_COLUMNS,
    SOURCE_COLUMNS,
    TARGET_COLUMNS,
    background_scores,
    empirical_weighted_mean,
    expression_fraction,
    filter_lr,
    normalize_gene,
    pick_column,
    prepare_magma,
    read_table,
    score_series,
    split_genes,
    standardize_lr,
    weighted_gene_set,
)


def enrich_tissue(
    lr_results: pd.DataFrame | str | Path,
    magma: pd.DataFrame | str | Path,
    lr_resource: pd.DataFrame | str | Path | None = None,
    magma_gene_col: str = "SYMBOL",
    magma_z_col: str = "ZSTAT",
    ligand_col: str | None = None,
    receptor_col: str | None = None,
    strength_col: str | None = None,
    score_mode: str = "abs_z",
    n_permutations: int = 5000,
    random_state: int | None = 0,
    replace: bool = False,
    dataset: str = "tissue",
) -> pd.DataFrame:
    """Test all filtered LR genes together at the tissue/dataset level."""

    lr = standardize_lr(lr_results, ligand_col=ligand_col, receptor_col=receptor_col)
    magma_table = prepare_magma(magma, gene_symbol_col=magma_gene_col, z_score_col=magma_z_col)
    scores = score_series(magma_table, score_mode=score_mode)
    genes, weights = weighted_gene_set(lr, range(len(lr)), strength_col=strength_col)
    foreground = scores.reindex(genes).dropna()
    keep = np.asarray([gene in set(foreground.index) for gene in genes])
    observed, null_mean, null_sd, enrichment_z, empirical_p = empirical_weighted_mean(
        foreground.to_numpy(dtype=float), weights[keep], background_scores(magma_table, lr_resource, score_mode), n_permutations, np.random.default_rng(random_state), replace=replace
    )
    return pd.DataFrame([{
        "dataset": dataset,
        "n_genes": int(len(foreground)),
        "observed_stat": observed,
        "null_mean": null_mean,
        "null_sd": null_sd,
        "enrichment_z": enrichment_z,
        "empirical_p": empirical_p,
    }])


def enrich_cell_type_pairs(
    lr_results: pd.DataFrame | str | Path,
    magma: pd.DataFrame | str | Path,
    lr_resource: pd.DataFrame | str | Path | None = None,
    magma_gene_col: str = "SYMBOL",
    magma_z_col: str = "ZSTAT",
    source_col: str | None = None,
    target_col: str | None = None,
    ligand_col: str | None = None,
    receptor_col: str | None = None,
    strength_col: str | None = None,
    score_mode: str = "abs_z",
    n_permutations: int = 5000,
    random_state: int | None = 0,
    replace: bool = False,
) -> pd.DataFrame:
    """Test LR-involved genes independently for each source-target pair."""

    lr = standardize_lr(lr_results, ligand_col=ligand_col, receptor_col=receptor_col, source_col=source_col, target_col=target_col)
    if not {"source", "target"}.issubset(lr.columns):
        raise KeyError("Cell-type pair enrichment requires source and target columns")
    magma_table = prepare_magma(magma, gene_symbol_col=magma_gene_col, z_score_col=magma_z_col)
    scores = score_series(magma_table, score_mode=score_mode)
    background = background_scores(magma_table, lr_resource, score_mode)
    rng = np.random.default_rng(random_state)
    rows = []
    for (source, target), group in lr.groupby(["source", "target"], sort=True, dropna=False):
        group = group.reset_index(drop=True)
        genes, weights = weighted_gene_set(group, range(len(group)), strength_col=strength_col)
        foreground = scores.reindex(genes).dropna()
        keep = np.asarray([gene in set(foreground.index) for gene in genes])
        if foreground.empty:
            continue
        observed, null_mean, null_sd, enrichment_z, empirical_p = empirical_weighted_mean(foreground.to_numpy(), weights[keep], background, n_permutations, rng, replace=replace)
        rows.append({
            "source": source,
            "target": target,
            "n_genes": int(len(foreground)),
            "observed_stat": observed,
            "null_mean": null_mean,
            "null_sd": null_sd,
            "enrichment_z": enrichment_z,
            "empirical_p": empirical_p,
        })
    if not rows:
        raise ValueError("No cell-type pair has genes matched to both LR and MAGMA")
    return pd.DataFrame(rows).sort_values(["empirical_p", "enrichment_z"], ascending=[True, False]).reset_index(drop=True)


def _pair_signal(ligand_z: float, receptor_z: float, mode: str) -> float:
    mode = mode.lower()
    if mode in {"geom_mean", "geometric_mean"}:
        return float(np.sqrt(max(ligand_z, 0.0) * max(receptor_z, 0.0)))
    if mode == "sum":
        return float(max(ligand_z, 0.0) + max(receptor_z, 0.0))
    if mode in {"prod", "product"}:
        return float(max(ligand_z, 0.0) * max(receptor_z, 0.0))
    raise ValueError("magma_mode must be geom_mean, sum, or prod")


def enrich_lr(
    lr_results: pd.DataFrame | str | Path,
    magma: pd.DataFrame | str | Path,
    lr_resource: pd.DataFrame | str | Path | None = None,
    magma_gene_col: str = "SYMBOL",
    magma_z_col: str = "ZSTAT",
    ligand_col: str | None = None,
    receptor_col: str | None = None,
    magma_mode: str = "geom_mean",
    n_permutations: int = 5000,
    random_state: int | None = 0,
) -> pd.DataFrame:
    """Test individual ligand-receptor pairs with empirical p-values."""

    lr = standardize_lr(lr_results, ligand_col=ligand_col, receptor_col=receptor_col)
    magma_table = prepare_magma(magma, gene_symbol_col=magma_gene_col, z_score_col=magma_z_col)
    z = magma_table.set_index("gene")["zscore"].astype(float)
    background = background_scores(magma_table, lr_resource, "z")
    rng = np.random.default_rng(random_state)
    cache = {}
    rows = []
    for _, row in lr.iterrows():
        ligand_genes = [gene for gene in split_genes(row["ligand"]) if gene in z.index]
        receptor_genes = [gene for gene in split_genes(row["receptor"]) if gene in z.index]
        if not ligand_genes or not receptor_genes:
            continue
        ligand_gene = max(ligand_genes, key=lambda gene: z.loc[gene])
        receptor_gene = max(receptor_genes, key=lambda gene: z.loc[gene])
        key = (ligand_gene, receptor_gene)
        if key not in cache:
            observed = _pair_signal(float(z.loc[ligand_gene]), float(z.loc[receptor_gene]), magma_mode)
            null = np.empty(n_permutations, dtype=float)
            for index in range(n_permutations):
                sampled = rng.choice(background, size=2, replace=len(background) < 2)
                null[index] = _pair_signal(float(sampled[0]), float(sampled[1]), magma_mode)
            null_mean = float(null.mean())
            null_sd = float(null.std(ddof=1)) if len(null) > 1 else 0.0
            cache[key] = (observed, (observed - null_mean) / null_sd if null_sd > 0 else np.nan, float((1 + np.count_nonzero(null >= observed)) / (n_permutations + 1)))
        observed, enrichment_z, empirical_p = cache[key]
        output = {
            "ligand": row["ligand"],
            "receptor": row["receptor"],
            "ligand_magma_gene": ligand_gene,
            "receptor_magma_gene": receptor_gene,
            "z_ligand": float(z.loc[ligand_gene]),
            "z_receptor": float(z.loc[receptor_gene]),
            "observed_stat": observed,
            "enrichment_z": enrichment_z,
            "empirical_p": empirical_p,
        }
        for column in ("source", "target"):
            if column in row.index:
                output[column] = row[column]
        rows.append(output)
    if not rows:
        raise ValueError("No LR rows matched ligand/receptor genes to MAGMA")
    result = pd.DataFrame(rows).sort_values(["empirical_p", "observed_stat"], ascending=[True, False]).reset_index(drop=True)
    columns = ["ligand", "receptor", "source", "target", "ligand_magma_gene", "receptor_magma_gene", "z_ligand", "z_receptor", "observed_stat", "enrichment_z", "empirical_p"]
    return result[[column for column in columns if column in result.columns]]


def _nichenet_long_table(nichenet: pd.DataFrame | str | Path, ligand_col: str | None = None, target_col: str | None = None, score_col: str | None = None) -> pd.DataFrame:
    table = read_table(nichenet)
    if ligand_col is not None and target_col is not None and score_col is not None:
        return pd.DataFrame({"ligand": table[ligand_col].map(normalize_gene), "target_gene": table[target_col].map(normalize_gene), "regulatory_potential": pd.to_numeric(table[score_col], errors="coerce")}).dropna()
    target = pick_column(table.columns, target_col, ("Target", "target", "target_gene", "gene"), "NicheNet target column")
    long_names = {"ligand", "source", "target_gene", "regulatory_potential", "weight", "score"}
    if ligand_col is not None or score_col is not None or len(set(table.columns).intersection(long_names)) >= 3:
        ligand = pick_column(table.columns, ligand_col, ("ligand", "source", "from"), "NicheNet ligand column")
        score = pick_column(table.columns, score_col, ("regulatory_potential", "weight", "score", "value"), "NicheNet score column")
        return pd.DataFrame({"ligand": table[ligand].map(normalize_gene), "target_gene": table[target].map(normalize_gene), "regulatory_potential": pd.to_numeric(table[score], errors="coerce")}).dropna()
    rows = []
    for column in table.columns:
        if column == target:
            continue
        ligand = normalize_gene(column)
        if ligand:
            rows.append(pd.DataFrame({"ligand": ligand, "target_gene": table[target].map(normalize_gene), "regulatory_potential": pd.to_numeric(table[column], errors="coerce")}))
    if not rows:
        raise ValueError("No ligand columns found in the NicheNet matrix")
    return pd.concat(rows, ignore_index=True).dropna()


def enrich_downstream(
    lr_results: pd.DataFrame | str | Path,
    adata: object,
    magma: pd.DataFrame | str | Path,
    nichenet: pd.DataFrame | str | Path,
    cell_type_col: str,
    magma_gene_col: str = "SYMBOL",
    magma_z_col: str = "ZSTAT",
    cellphone_p_cutoff: float | None = 0.05,
    lr_filter_query: str | None = None,
    top_n_targets: int = 300,
    expression_fraction_cutoff: float = 0.05,
    min_genes: int = 10,
    n_permutations: int = 5000,
    random_state: int | None = 0,
    ligand_col: str | None = None,
    receptor_col: str | None = None,
    source_col: str | None = None,
    target_col: str | None = None,
    nichenet_ligand_col: str | None = None,
    nichenet_target_col: str | None = None,
    nichenet_score_col: str | None = None,
) -> pd.DataFrame:
    """Test expressed NicheNet targets for each selected LR context.

    The default LR filter is ``cellphone_pvals < 0.05``. Use
    ``lr_filter_query`` for another LIANA/Liana+ criterion.
    """

    if lr_filter_query is not None:
        selected = filter_lr(lr_results, filter_query=lr_filter_query)
    elif cellphone_p_cutoff is not None:
        raw_lr = read_table(lr_results)
        p_column = pick_column(raw_lr.columns, None, CELL_PHONE_COLUMNS, "CellPhoneDB p-value column")
        selected = filter_lr(raw_lr, filter_col=p_column, threshold=cellphone_p_cutoff)
    else:
        selected = read_table(lr_results)
    lr = standardize_lr(selected, ligand_col=ligand_col, receptor_col=receptor_col, source_col=source_col, target_col=target_col)
    if not {"source", "target"}.issubset(lr.columns):
        raise KeyError("Downstream enrichment requires source and target columns")
    if "cellphone_pvals" not in lr.columns:
        for column in CELL_PHONE_COLUMNS:
            if column in lr.columns:
                lr["cellphone_pvals"] = pd.to_numeric(lr[column], errors="coerce")
                break
    network = _nichenet_long_table(nichenet, ligand_col=nichenet_ligand_col, target_col=nichenet_target_col, score_col=nichenet_score_col)
    magma_table = prepare_magma(magma, gene_symbol_col=magma_gene_col, z_score_col=magma_z_col)
    scores = score_series(magma_table, score_mode="abs_z")
    background = scores.reindex(network["target_gene"].dropna().unique()).dropna().to_numpy(dtype=float)
    if background.size == 0:
        raise ValueError("No NicheNet target genes overlap MAGMA")
    rng = np.random.default_rng(random_state)
    rows = []
    for _, row in lr.iterrows():
        targets = network[network["ligand"].isin(split_genes(row["ligand"]))].copy()
        if targets.empty:
            continue
        targets = targets.groupby("target_gene", as_index=False)["regulatory_potential"].max().sort_values("regulatory_potential", ascending=False).head(top_n_targets)
        fractions = expression_fraction(adata, targets["target_gene"], cell_type_col=cell_type_col, cell_types=[row["target"]])
        if not fractions.empty:
            targets = targets.merge(fractions.rename(columns={"cell_type": "target"}), left_on="target_gene", right_on="gene", how="left")
            targets = targets[targets["expr_frac"].fillna(0) >= expression_fraction_cutoff]
        targets = targets[targets["target_gene"].isin(scores.index)].drop_duplicates("target_gene")
        if len(targets) < min_genes:
            continue
        genes = targets["target_gene"].tolist()
        weights = targets["regulatory_potential"].to_numpy(dtype=float)
        observed, null_mean, null_sd, enrichment_z, empirical_p = empirical_weighted_mean(scores.loc[genes].to_numpy(), weights, background, n_permutations, rng)
        output = {
            "source": row["source"],
            "target": row["target"],
            "ligand": row["ligand"],
            "receptor": row["receptor"],
            "n_genes": len(genes),
            "observed_stat": observed,
            "null_mean": null_mean,
            "null_sd": null_sd,
            "enrichment_z": enrichment_z,
            "empirical_p": empirical_p,
        }
        if "cellphone_pvals" in row.index:
            output["cellphone_pvals"] = row["cellphone_pvals"]
        rows.append(output)
    if not rows:
        raise ValueError("No downstream context had enough expressed NicheNet targets matched to MAGMA")
    return pd.DataFrame(rows).sort_values(["empirical_p", "enrichment_z"], ascending=[True, False]).reset_index(drop=True)


__all__ = ["enrich_tissue", "enrich_cell_type_pairs", "enrich_lr", "enrich_downstream"]
