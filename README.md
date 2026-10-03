# molplumber

A small command-line tool and Python function for checking **train/test leakage** in molecular datasets.

For every test molecule, `molplumber` reports:

- the **maximum Tanimoto similarity** to any training molecule, and which training molecule that is;
- whether the molecule is an **exact match** (canonical SMILES) to a training molecule;
- its **Bemis-Murcko scaffold**, and whether that scaffold also occurs in the training set.

Reported model performance depends on how close the test set is to the training set. This tool does not fix a split. It measures how close the two sets are, so you can state it alongside your metrics.

## Installation

Inside your project directory where `pyproject.toml` is located:

### Using `uv`

```bash
uv pip install -e .

```

### Using standard `pip`

```bash
pip install -e .

```

Requires Python ≥ 3.9, RDKit, pandas and numpy.

## Usage

### Command line

```bash
molplumber train.csv test.csv --smiles-col smiles -o leakage_report.csv
```

| Option | Meaning |
|---|---|
| `--smiles-col` | SMILES column name (default `smiles`); override per file with `--train-smiles-col` / `--test-smiles-col` |
| `--radius`, `--n-bits` | Morgan fingerprint settings (default 2, 2048) |
| `--generic-scaffold` | Compare generic scaffolds (atom and bond types ignored) instead of standard ones |
| `--thresholds` | Similarity cutoffs used in the printed summary (default `0.4 0.6 0.8 1.0`) |
| `--summary-json` | Also save the summary as JSON |

Files ending in `.tsv` or `.txt` are read as tab-separated; everything else as comma-separated.

### Python

```python
import pandas as pd
from molplumber import assess_leakage, summarize

train = pd.read_csv("train.csv")
test = pd.read_csv("test.csv")

report = assess_leakage(train["smiles"], test["smiles"])
print(summarize(report))
```

## Output

`leakage_report.csv` has one row per test molecule, in input order.

| Column | Description |
|---|---|
| `test_index` | Row position in the test file (0-based) |
| `smiles` | Input SMILES |
| `valid` | `False` if RDKit could not parse the SMILES; other fields are then empty |
| `max_tanimoto` | Highest Tanimoto similarity to any training molecule |
| `nearest_train_index`, `nearest_train_smiles` | The training molecule that gave `max_tanimoto` |
| `exact_match` | Canonical SMILES (including stereochemistry) also present in the training set |
| `scaffold` | Murcko scaffold SMILES (stereochemistry removed); empty for acyclic molecules |
| `scaffold_in_train` | Scaffold occurs in the training set; empty for acyclic molecules |
| `n_train_with_scaffold` | Number of training molecules with that scaffold |

The printed summary gives the distribution of `max_tanimoto`, the fraction of test molecules at or above each threshold, the exact-match fraction, and the fraction of test scaffolds seen in training.

## Interpreting the results

- A test set whose molecules mostly have close training analogs will tend to give optimistic performance estimates for new chemotypes. A low-similarity, low-overlap test set is a harder and often more realistic evaluation.
- The default thresholds are arbitrary reporting cutoffs, not established pass/fail limits. Which values are meaningful depends on your fingerprint and your task, so report the full distribution where you can.
- Tanimoto values are specific to the fingerprint. Changing `--radius` or `--n-bits` changes them, so state the settings used.

## Limitations

- **Acyclic molecules** have no Murcko scaffold. They are excluded from the scaffold-overlap fraction and counted separately as `n_acyclic_test`, since they would otherwise all share an empty scaffold.
- **Exact match** uses canonical SMILES from a single RDKit parse. Different tautomers, charge states, or salts of the same compound count as different molecules. Standardize your data first if that matters.
- **Scaling.** Similarity is computed by brute force (every test molecule against every training molecule). This is fine for tens of thousands of molecules; it has not been benchmarked on much larger sets.
- Invalid SMILES are reported and skipped, not repaired.

## Using the report in R

```r
library(ggplot2)

rep <- read.csv("leakage_report.csv")
ggplot(subset(rep, valid), aes(max_tanimoto)) +
  geom_histogram(bins = 30) +
  labs(x = "Max Tanimoto similarity to training set", y = "Test molecules")
```

## Tests

```bash
pip install -e ".[dev]"
pytest
```

---

## License

MIT

## Disclaimer

Google Gemini 3.8 Flash was used to tidy up the code and generate the README.md file.
