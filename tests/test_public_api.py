import pandas as pd

import commgwas


def example_inputs():
    magma = pd.DataFrame(
        {
            "SYMBOL": ["LIG1", "REC1", "LIG2", "REC2", "GENE1", "GENE2", "GENE3"],
            "ZSTAT": [3.0, 2.5, 2.0, 1.5, 0.3, 0.2, 0.1],
        }
    )
    lr = pd.DataFrame(
        {
            "ligand": ["LIG1", "LIG2"],
            "receptor": ["REC1", "REC2"],
            "source": ["A", "B"],
            "target": ["B", "A"],
            "cellphone_pvals": [0.01, 0.02],
        }
    )
    resource = pd.DataFrame(
        {
            "ligand": ["LIG1", "LIG2", "GENE1"],
            "receptor": ["REC1", "REC2", "GENE2"],
        }
    )
    return magma, lr, resource


def test_public_enrichment_outputs_do_not_create_fdr():
    magma, lr, resource = example_inputs()
    filtered = commgwas.filter_lr(lr)
    for result in (
        commgwas.enrich_tissue(filtered, magma, lr_resource=resource, n_permutations=20),
        commgwas.enrich_cell_type_pairs(filtered, magma, lr_resource=resource, n_permutations=20),
        commgwas.enrich_lr(filtered, magma, lr_resource=resource, n_permutations=20),
    ):
        assert isinstance(result, pd.DataFrame)
        assert not any("fdr" in str(column).lower() for column in result.columns)
        assert "empirical_p" in result.columns
