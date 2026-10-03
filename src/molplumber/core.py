"""Core functions for assessing train/test leakage in molecular datasets.

For every test molecule the report gives:
  * the maximum Tanimoto similarity to any training molecule,
  * the nearest training molecule,
  * whether the test molecule is an exact match (canonical SMILES) to a training molecule,
  * its Bemis-Murcko scaffold and whether that scaffold occurs in the training set.
"""
from __future__ import annotations

from collections import Counter
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

DEFAULT_THRESHOLDS = (0.4, 0.6, 0.8, 1.0)


def _parse(smiles) -> Chem.Mol | None:
    """Return an RDKit Mol, or None if the input is not a parseable SMILES string."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    return Chem.MolFromSmiles(smiles)


def murcko_scaffold(mol: Chem.Mol, generic: bool = False, chirality: bool = False) -> str | None:
    """Bemis-Murcko scaffold as SMILES.

    Returns None for acyclic molecules (empty scaffold) or if scaffold
    generation fails. Acyclic molecules are excluded from scaffold overlap
    statistics, because they would all share the same empty scaffold.
    """
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        if generic:
            scaf = MurckoScaffold.MakeScaffoldGeneric(scaf)
        smi = Chem.MolToSmiles(scaf, isomericSmiles=chirality)
    except Exception:
        return None
    return smi or None


def assess_leakage(
    train_smiles: Iterable[str],
    test_smiles: Iterable[str],
    radius: int = 2,
    n_bits: int = 2048,
    generic_scaffold: bool = False,
) -> pd.DataFrame:
    """Compare each test molecule against the training set.

    Parameters
    ----------
    train_smiles, test_smiles : iterables of SMILES strings.
    radius, n_bits : Morgan (ECFP-like) bit-vector fingerprint settings.
    generic_scaffold : if True, compare generic (atom/bond-type-agnostic)
        scaffolds instead of standard Bemis-Murcko scaffolds.

    Returns
    -------
    DataFrame with one row per test molecule, in input order.
    """
    train_smiles = list(train_smiles)
    test_smiles = list(test_smiles)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)

    # --- training set ---
    train_idx, train_fps, train_can, train_scaf = [], [], [], []
    for i, smi in enumerate(train_smiles):
        mol = _parse(smi)
        if mol is None:
            continue
        train_idx.append(i)
        train_fps.append(gen.GetFingerprint(mol))
        train_can.append(Chem.MolToSmiles(mol))
        train_scaf.append(murcko_scaffold(mol, generic=generic_scaffold))
    if not train_fps:
        raise ValueError("No valid molecules in the training set.")

    can_set = set(train_can)
    scaf_counts = Counter(s for s in train_scaf if s)

    # --- test set ---
    rows = []
    for j, smi in enumerate(test_smiles):
        mol = _parse(smi)
        if mol is None:
            rows.append({"test_index": j, "smiles": smi, "valid": False})
            continue
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(gen.GetFingerprint(mol), train_fps))
        k = int(sims.argmax())
        scaf = murcko_scaffold(mol, generic=generic_scaffold)
        rows.append(
            {
                "test_index": j,
                "smiles": smi,
                "valid": True,
                "max_tanimoto": float(sims[k]),
                "nearest_train_index": train_idx[k],
                "nearest_train_smiles": train_smiles[train_idx[k]],
                "exact_match": Chem.MolToSmiles(mol) in can_set,
                "scaffold": scaf,
                "scaffold_in_train": (scaf in scaf_counts) if scaf else None,
                "n_train_with_scaffold": scaf_counts.get(scaf, 0) if scaf else None,
            }
        )

    df = pd.DataFrame(rows)
    df["valid"] = df["valid"].astype(bool)
    for col in ("exact_match", "scaffold_in_train"):
        if col in df:
            df[col] = df[col].astype("boolean")
    for col in ("nearest_train_index", "n_train_with_scaffold"):
        if col in df:
            df[col] = df[col].astype("Int64")
    return df


def summarize(df: pd.DataFrame, thresholds: Sequence[float] = DEFAULT_THRESHOLDS) -> dict:
    """Aggregate a per-molecule report into summary statistics."""
    valid = df[df["valid"]]
    out: dict = {
        "n_test": int(len(df)),
        "n_invalid_test": int((~df["valid"]).sum()),
        "n_valid_test": int(len(valid)),
    }
    if valid.empty:
        return out

    sims = valid["max_tanimoto"]
    out["max_tanimoto_min"] = float(sims.min())
    out["max_tanimoto_q25"] = float(sims.quantile(0.25))
    out["max_tanimoto_median"] = float(sims.median())
    out["max_tanimoto_q75"] = float(sims.quantile(0.75))
    out["max_tanimoto_max"] = float(sims.max())
    out["frac_exact_match"] = float(valid["exact_match"].mean())
    for t in thresholds:
        out[f"frac_max_tanimoto_ge_{t}"] = float((sims >= t).mean())

    cyclic = valid[valid["scaffold"].notna()]
    out["n_acyclic_test"] = int(len(valid) - len(cyclic))
    out["frac_scaffold_in_train"] = (
        float(cyclic["scaffold_in_train"].mean()) if len(cyclic) else None
    )
    return out
