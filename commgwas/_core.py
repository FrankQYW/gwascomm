"""Internal utilities for the standalone commgwas package."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

LIGAND_COLUMNS = ("ligand", "Ligand", "ligand_complex", "source_genesymbol", "from")
RECEPTOR_COLUMNS = ("receptor", "Receptor", "receptor_complex", "target_genesymbol", "to")
SOURCE_COLUMNS = ("source", "sender", "Sender", "celltype_sender")
TARGET_COLUMNS = ("target", "receiver", "Receiver", "celltype_receiver")
CELL_PHONE_COLUMNS = ("cellphone_pvals", "cellphone_pval", "cellphone_p", "p_CCC", "pvalue", "p_value", "p")


def read_table(table: pd.DataFrame | str | Path) -> pd.DataFrame:
    if isinstance(table, pd.DataFrame):
        return table.copy()
    path = Path(table)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t")
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def pick_column(columns: Iterable[object], requested: str | None, candidates: Sequence[str], label: str) -> str:
    names = [str(column) for column in columns]
    if requested is not None:
        if requested not in names:
            raise KeyError(f"{label} `{requested}` was not found")
        return requested
    for candidate in candidates:
        if candidate in names:
            return candidate
    lower = {name.lower(): name for name in names}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    raise KeyError(f"Could not find {label}; tried {', '.join(candidates)}")


def normalize_gene(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    return text.upper() or None


def split_genes(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    for separator in ("^", "&", ";", "|", ","):
        text = text.replace(separator, "_")
    return [gene for gene in (normalize_gene(part) for part in text.split("_")) if gene]


def prepare_magma(
    magma: pd.DataFrame | str | Path,
    gene_symbol_col: str = "SYMBOL",
    z_score_col: str = "ZSTAT",
    p_value_col: str | None = None,
) -> pd.DataFrame:
    """Normalize a MAGMA table to only the columns needed by commgwas."""

    table = read_table(magma)
    for column, label in ((gene_symbol_col, "gene symbol"), (z_score_col, "z-score")):
        if column not in table.columns:
            raise KeyError(f"MAGMA {label} column `{column}` was not found")
    out = pd.DataFrame(
        {
            "gene": table[gene_symbol_col].map(normalize_gene),
            "zscore": pd.to_numeric(table[z_score_col], errors="coerce"),
        }
    )
    if p_value_col is not None:
        if p_value_col not in table.columns:
            raise KeyError(f"MAGMA p-value column `{p_value_col}` was not found")
        out["pvalue"] = pd.to_numeric(table[p_value_col], errors="coerce")
    out = out.dropna(subset=["gene", "zscore"])
    out = out[np.isfinite(out["zscore"].to_numpy(dtype=float))]
    if out.empty:
        raise ValueError("MAGMA contains no valid gene/z-score rows")
    return out.sort_values("zscore", ascending=False).drop_duplicates("gene").reset_index(drop=True)


def standardize_lr(
    lr_results: pd.DataFrame | str | Path,
    ligand_col: str | None = None,
    receptor_col: str | None = None,
    source_col: str | None = None,
    target_col: str | None = None,
) -> pd.DataFrame:
    table = read_table(lr_results)
    if table.empty:
        raise ValueError("The LR result table is empty")
    ligand = pick_column(table.columns, ligand_col, LIGAND_COLUMNS, "ligand column")
    receptor = pick_column(table.columns, receptor_col, RECEPTOR_COLUMNS, "receptor column")
    out = table.copy()
    out["ligand"] = out[ligand].where(out[ligand].notna(), np.nan).astype(object)
    out["receptor"] = out[receptor].where(out[receptor].notna(), np.nan).astype(object)
    if source_col is not None or any(column in out.columns for column in SOURCE_COLUMNS):
        source = pick_column(out.columns, source_col, SOURCE_COLUMNS, "source cell type column")
        out["source"] = out[source].astype(str)
    if target_col is not None or any(column in out.columns for column in TARGET_COLUMNS):
        target = pick_column(out.columns, target_col, TARGET_COLUMNS, "target cell type column")
        out["target"] = out[target].astype(str)
    out = out.dropna(subset=["ligand", "receptor"]).reset_index(drop=True)
    if out.empty:
        raise ValueError("No usable ligand-receptor rows remain")
    return out


def filter_lr(
    lr_results: pd.DataFrame | str | Path,
    filter_query: str | None = None,
    filter_func: Callable[[pd.DataFrame], Iterable[bool]] | None = None,
    filter_col: str | None = "cellphone_pvals",
    filter_op: str = "<",
    threshold: float = 0.05,
) -> pd.DataFrame:
    """Filter LIANA rows; the default is ``cellphone_pvals < 0.05``."""

    table = read_table(lr_results)
    if filter_query is not None:
        result = table.query(filter_query, engine="python")
    elif filter_func is not None:
        result = table.loc[np.asarray(list(filter_func(table)), dtype=bool)]
    elif filter_col is None:
        result = table
    else:
        if filter_col not in table.columns:
            raise KeyError(f"LR filter column `{filter_col}` was not found")
        values = pd.to_numeric(table[filter_col], errors="coerce")
        operators = {"<": values < threshold, "<=": values <= threshold, ">": values > threshold, ">=": values >= threshold, "==": values == threshold}
        if filter_op not in operators:
            raise ValueError("filter_op must be one of <, <=, >, >=, ==")
        result = table.loc[operators[filter_op]]
    if result.empty:
        raise ValueError("No LR rows remain after filtering")
    return result.copy().reset_index(drop=True)


def load_lr_resource(resource: pd.DataFrame | str | Path | None = None) -> pd.DataFrame:
    if resource is None or str(resource).lower() in {"consensus", "liana", "default"}:
        try:
            import liana as li
        except ImportError as exc:
            raise ImportError("Install `liana` to use the default consensus LR resource") from exc
        table = li.resource.select_resource("consensus")
    else:
        table = read_table(resource)
    ligand = pick_column(table.columns, None, LIGAND_COLUMNS, "LR resource ligand column")
    receptor = pick_column(table.columns, None, RECEPTOR_COLUMNS, "LR resource receptor column")
    return pd.DataFrame({"ligand": table[ligand].astype(str), "receptor": table[receptor].astype(str)}).drop_duplicates().reset_index(drop=True)


def score_series(magma: pd.DataFrame, score_mode: str = "abs_z") -> pd.Series:
    z = magma.set_index("gene")["zscore"].astype(float)
    mode = score_mode.lower()
    if mode in {"abs_z", "abs", "absolute_z"}:
        values = z.abs()
    elif mode in {"z", "raw_z"}:
        values = z
    elif mode in {"positive_z", "positive"}:
        values = z.clip(lower=0)
    else:
        raise ValueError("score_mode must be abs_z, z, or positive_z")
    values.name = mode
    return values.replace([np.inf, -np.inf], np.nan).dropna()


def pair_strength(table: pd.DataFrame, strength_col: str | None = None) -> np.ndarray:
    if strength_col is not None:
        if strength_col not in table.columns:
            raise KeyError(f"Strength column `{strength_col}` was not found")
        return np.abs(pd.to_numeric(table[strength_col], errors="coerce").fillna(0).to_numpy(dtype=float))
    if {"lr_means", "cellphone_pvals"}.issubset(table.columns):
        magnitude = pd.to_numeric(table["lr_means"], errors="coerce").fillna(0).clip(lower=0)
        maximum = float(magnitude.max())
        if maximum > 0:
            magnitude = magnitude / maximum
        pvalues = pd.to_numeric(table["cellphone_pvals"], errors="coerce").clip(lower=1e-300, upper=1).fillna(1)
        return (magnitude * -np.log10(pvalues)).to_numpy(dtype=float)
    return np.ones(len(table), dtype=float)


def weighted_gene_set(table: pd.DataFrame, indices: Iterable[int], strength_col: str | None = None) -> tuple[list[str], np.ndarray]:
    strengths = pair_strength(table, strength_col=strength_col)
    weights_by_gene: dict[str, float] = {}
    for index in indices:
        genes = set(split_genes(table.iloc[index]["ligand"]) + split_genes(table.iloc[index]["receptor"]))
        for gene in genes:
            weights_by_gene[gene] = weights_by_gene.get(gene, 0.0) + max(float(strengths[index]), 1e-12)
    genes = sorted(weights_by_gene)
    return genes, np.asarray([weights_by_gene[gene] for gene in genes], dtype=float)


def empirical_weighted_mean(
    observed_scores: np.ndarray,
    weights: np.ndarray,
    background_scores: np.ndarray,
    n_permutations: int,
    rng: np.random.Generator,
    replace: bool = False,
) -> tuple[float, float, float, float, float]:
    if len(observed_scores) == 0:
        raise ValueError("The foreground gene set is empty")
    if not replace and len(observed_scores) > len(background_scores):
        raise ValueError("Foreground is larger than the background; use replace=True or provide a larger universe")
    observed = float(np.average(observed_scores, weights=weights))
    null = np.empty(n_permutations, dtype=float)
    for index in range(n_permutations):
        sampled = rng.choice(background_scores, size=len(weights), replace=replace)
        null[index] = np.average(sampled, weights=weights)
    null_mean = float(null.mean())
    null_sd = float(null.std(ddof=1)) if len(null) > 1 else 0.0
    empirical_p = float((1 + np.count_nonzero(null >= observed)) / (n_permutations + 1))
    enrichment_z = (observed - null_mean) / null_sd if null_sd > 0 else np.nan
    return observed, null_mean, null_sd, enrichment_z, empirical_p


def background_scores(magma: pd.DataFrame, resource: pd.DataFrame | str | Path | None, score_mode: str) -> np.ndarray:
    scores = score_series(magma, score_mode=score_mode)
    if resource is None:
        return scores.to_numpy(dtype=float)
    lr = load_lr_resource(resource)
    genes = sorted({gene for value in pd.concat([lr["ligand"], lr["receptor"]]) for gene in split_genes(value)})
    values = scores.reindex(genes).dropna()
    if values.empty:
        raise ValueError("No LR-universe genes overlap the MAGMA table")
    return values.to_numpy(dtype=float)


def expression_fraction(adata: object, genes: Iterable[str], cell_type_col: str, cell_types: Iterable[str] | None = None, feature_col: str = "feature_name", layer: str | None = None) -> pd.DataFrame:
    try:
        import anndata as ad
        from scipy import sparse
    except ImportError as exc:
        raise ImportError("Install `anndata` and `scipy` for expression filtering") from exc
    close = False
    data = adata
    if isinstance(adata, (str, Path)):
        data = ad.read_h5ad(adata, backed="r")
        close = True
    try:
        if cell_type_col not in data.obs.columns:
            raise KeyError(f"Cell type column `{cell_type_col}` was not found in AnnData.obs")
        wanted = sorted({gene for gene in (normalize_gene(value) for value in genes) if gene})
        symbols = data.var[feature_col].astype(str).map(normalize_gene).to_numpy() if feature_col in data.var.columns else np.asarray([normalize_gene(value) for value in data.var_names])
        positions = {gene: index for index, gene in enumerate(symbols) if gene in wanted}
        values = data.obs[cell_type_col].astype(str)
        selected_types = sorted(set(values)) if cell_types is None else [str(value) for value in cell_types]
        rows = []
        for cell_type in selected_types:
            cell_indices = np.flatnonzero(values.to_numpy() == cell_type)
            if cell_indices.size == 0:
                continue
            selected_genes = [gene for gene in wanted if gene in positions]
            matrix = data[cell_indices, [positions[gene] for gene in selected_genes]].layers[layer] if layer is not None else data[cell_indices, [positions[gene] for gene in selected_genes]].X
            fractions = np.asarray((matrix > 0).mean(axis=0)).ravel() if sparse.issparse(matrix) else (np.asarray(matrix) > 0).mean(axis=0)
            rows.extend({"cell_type": cell_type, "gene": gene, "expr_frac": float(value)} for gene, value in zip(selected_genes, fractions))
        return pd.DataFrame(rows)
    finally:
        if close and getattr(data, "file", None) is not None:
            data.file.close()
