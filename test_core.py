import math

import pandas as pd

from molplumber import assess_leakage, summarize

TRAIN = ["c1ccccc1", "CCO", "c1ccc2ccccc2c1"]
TEST = ["c1ccccc1", "Cc1ccccc1", "CCCO", "not_a_smiles", "C1CCNCC1"]


def _report():
    return assess_leakage(TRAIN, TEST)


def test_exact_match_and_similarity():
    r = _report().set_index("test_index")
    assert r.loc[0, "exact_match"]
    assert math.isclose(r.loc[0, "max_tanimoto"], 1.0)
    assert not r.loc[1, "exact_match"]
    assert r.loc[1, "max_tanimoto"] < 1.0


def test_scaffold_overlap():
    r = _report().set_index("test_index")
    assert r.loc[1, "scaffold"] == "c1ccccc1"
    assert r.loc[1, "scaffold_in_train"]       # toluene -> benzene, present in train
    assert not r.loc[4, "scaffold_in_train"]   # piperidine, absent from train


def test_acyclic_has_no_scaffold():
    r = _report().set_index("test_index")
    assert pd.isna(r.loc[2, "scaffold"])
    assert pd.isna(r.loc[2, "scaffold_in_train"])


def test_invalid_smiles_flagged():
    r = _report().set_index("test_index")
    assert not r.loc[3, "valid"]


def test_summary_counts():
    s = summarize(_report())
    assert s["n_test"] == 5
    assert s["n_invalid_test"] == 1
    assert s["n_acyclic_test"] == 1

