# commgwas

`commgwas` is a standalone, Scanpy-style Python package for testing whether
cell-cell communication programs are enriched for MAGMA gene-level GWAS
signals.

The package is self-contained under `CommGWAS/package/`. Its public functions accept
pandas DataFrames, table paths, and AnnData objects, and return pandas
DataFrames. It has no console entry point, command dispatcher, or command-line
workflow.

## Public API

```python
import commgwas

magma = commgwas.prepare_magma(
    magma_df,
    gene_symbol_col="SYMBOL",
    z_score_col="ZSTAT",
)

lr_raw = commgwas.infer_lr(adata, cell_type_col="predicted.celltype.l2")
lr_sig = commgwas.filter_lr(lr_raw)  # cellphone_pvals < 0.05

tissue = commgwas.enrich_tissue(lr_sig, magma_df)
cell_pairs = commgwas.enrich_cell_type_pairs(lr_sig, magma_df)
lr_level = commgwas.enrich_lr(lr_sig, magma_df)
```

The four enrichment functions return only the requested analysis columns and
raw empirical p-values. No `fdr`, `qvalue`, or Bonferroni column is created.
Users can apply their preferred multiple-testing procedure externally to the
returned `empirical_p` column.

### Filtering

`filter_lr()` defaults to the CellPhoneDB criterion:

```python
lr_sig = commgwas.filter_lr(lr_raw, threshold=0.05)
```

Other LIANA or LIANA+ criteria can be supplied without changing the package:

```python
lr_sig = commgwas.filter_lr(
    lr_raw,
    filter_query="cellphone_pvals < 0.01 and aggregate_rank < 0.2",
)
```

Tissue and cell-type-pair enrichment assumes that the LR table has already
been filtered. Downstream enrichment can either apply the default CellPhoneDB
filter or accept an explicit `lr_filter_query`.

### Downstream targets

```python
nichenet = "../resource/ligand_target_matrix_nsga2r_final.tsv"

downstream = commgwas.enrich_downstream(
    lr_results=lr_sig,
    adata=adata,
    magma=magma_df,
    nichenet=nichenet,
    cell_type_col="predicted.celltype.l2",
    cellphone_p_cutoff=None,  # lr_sig was filtered explicitly above
    min_genes=10,
)
```

The NicheNet input may be a wide matrix with a `Target` column or a long table
with ligand, target, and regulatory-potential columns.

## Example data

The tutorial `tutorials/commgwas_blood_sle.ipynb` uses:

- MAGMA: `magma/magma_cao/DisCP0423_1.magma.genes.results.csv`
- AnnData: `data/Blood_Hao/qc.h5ad`
- NicheNet: `resource/ligand_target_matrix_nsga2r_final.tsv`
- Cell-type metadata column: `predicted.celltype.l2`

The tutorial covers tissue-level, cell-type-pair, LR-pair, and NicheNet
 downstream enrichment results.

## Installation

Install the package from this directory with the base dependencies. Install
the optional `singlecell` dependencies when using `infer_lr()` or AnnData-based
expression filtering.

```python
# In a notebook environment:
%pip install -e ".[singlecell]"
```
