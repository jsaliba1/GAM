import pandas as pd

from gam.models.group_association import GroupAssociationModel


def _build_model() -> GroupAssociationModel:
    sample_ids = [f"U{i}" for i in range(1, 21)]

    bug_df = pd.DataFrame({"ID": sample_ids})
    drug_df = pd.DataFrame(
        {
            "UNIQUEID": sample_ids,
            "DRUG_A": ["S"] * 10 + ["R"] * 10,
            "DRUG_B": ["S"] * 20,
        }
    ).set_index("UNIQUEID")

    mutation_rows = []
    for uid in sample_ids:
        mutation_rows.append(
            {"UNIQUEID": uid, "MUTATION": "shared", "GENE": "g0", "IS_SYNONYMOUS": False}
        )
    for uid in sample_ids[10:]:
        mutation_rows.append(
            {"UNIQUEID": uid, "MUTATION": "mut_x", "GENE": "g1", "IS_SYNONYMOUS": False}
        )

    mutation_df = pd.DataFrame(mutation_rows)
    return GroupAssociationModel(mutations=mutation_df, drug_profiles=drug_df, bugs=bug_df)


def test_group_association_pipeline_returns_scores() -> None:
    model = _build_model()
    model.group_bugs(mask_probability=0.0, remove_n=0, seed=0)
    pvals = model.analyze_resistance()
    scored = model.score_resistance()

    assert pvals
    assert "DRUG_A" in scored.columns
