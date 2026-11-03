"""Scanpy-style cell-cell communication enrichment for MAGMA signals."""

from ._core import filter_lr, load_lr_resource, prepare_magma
from .enrichment import enrich_cell_type_pairs, enrich_downstream, enrich_lr, enrich_tissue
from .singlecell import infer_lr

__version__ = "0.1.0"

__all__ = [
    "prepare_magma",
    "load_lr_resource",
    "infer_lr",
    "filter_lr",
    "enrich_tissue",
    "enrich_cell_type_pairs",
    "enrich_lr",
    "enrich_downstream",
]
