"""Scanpy/AnnData-facing communication inference."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ._core import load_lr_resource


def infer_lr(
    adata: object,
    cell_type_col: str,
    lr_resource: pd.DataFrame | str | Path | None = None,
    min_cells: int = 20,
    min_expr_fraction: float = 0.05,
    n_permutations: int = 1000,
    use_raw: bool = False,
    layer: str | None = None,
    n_jobs: int = 1,
    random_state: int = 0,
) -> pd.DataFrame:
    """Run LIANA consensus aggregation on an AnnData object or h5ad path."""

    try:
        import anndata as ad
        import liana as li
    except ImportError as exc:
        raise ImportError("Install the `singlecell` extras to run infer_lr") from exc
    data = ad.read_h5ad(adata) if isinstance(adata, (str, Path)) else adata.copy()
    if cell_type_col not in data.obs.columns:
        raise KeyError(f"Cell type column `{cell_type_col}` was not found in AnnData.obs")
    kwargs = {
        "groupby": cell_type_col,
        "expr_prop": min_expr_fraction,
        "min_cells": min_cells,
        "key_added": "commgwas_liana",
        "use_raw": use_raw,
        "n_perms": n_permutations,
        "seed": random_state,
        "n_jobs": n_jobs,
        "inplace": True,
        "verbose": False,
    }
    if layer is not None:
        kwargs["layer"] = layer
    if lr_resource is None or str(lr_resource).lower() in {"consensus", "liana", "default"}:
        kwargs["resource_name"] = "consensus"
    else:
        kwargs["resource"] = load_lr_resource(lr_resource)
    li.mt.rank_aggregate(data, **kwargs)
    if "commgwas_liana" not in data.uns:
        raise ValueError("LIANA did not produce a result table")
    result = data.uns["commgwas_liana"].copy()
    if "ligand" not in result.columns and "ligand_complex" in result.columns:
        result["ligand"] = result["ligand_complex"]
    if "receptor" not in result.columns and "receptor_complex" in result.columns:
        result["receptor"] = result["receptor_complex"]
    if "interaction" not in result.columns and {"ligand", "receptor"}.issubset(result.columns):
        result["interaction"] = result["ligand"].astype(str) + "^" + result["receptor"].astype(str)
    return result.reset_index(drop=True)


__all__ = ["infer_lr"]
