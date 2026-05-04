import pandas as pd

from gam.config import GamSchema
from gam.models.group_association import GroupAssociationModel


def test_group_association_supports_custom_column_mapping() -> None:
    sample_ids = [f"S{i}" for i in range(1, 21)]
    mutations = []
    for sample_id in sample_ids:
        mutations.append(
            {"sample": sample_id, "variant": "shared", "gene_name": "g0", "syn": False}
        )
    for sample_id in sample_ids[10:]:
        mutations.append({"sample": sample_id, "variant": "mut_x", "gene_name": "g1", "syn": False})

    mutation_df = pd.DataFrame(mutations)
    bug_df = pd.DataFrame({"sample_name": sample_ids})
    drug_df = pd.DataFrame(
        {"sample_name": sample_ids, "drug_a": ["0"] * 10 + ["1"] * 10}
    ).set_index("sample_name")

    schema = GamSchema(
        mutation_sample_id_column="sample",
        mutation_name_column="variant",
        mutation_gene_column="gene_name",
        mutation_is_synonymous_column="syn",
        bug_sample_id_column="sample_name",
        drug_sample_id_column="sample_name",
        susceptible_label="0",
        resistant_label="1",
    )

    model = GroupAssociationModel(
        mutations=mutation_df,
        drug_profiles=drug_df,
        bugs=bug_df,
        schema=schema,
    )
    model.group_bugs(mask_probability=0.0, remove_n=0, seed=0)
    pvals = model.analyze_resistance()
    assert pvals
