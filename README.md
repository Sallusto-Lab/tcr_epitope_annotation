# tcr_epitope_annotation

Utilities for matching T-cell receptor (TCR) clonotypes against three
TCR–epitope databases:

- VDJdb
- IEDB
- McPAS-TCR

The code supports exact CDR3 matching and optional Levenshtein-distance
matching.

## Requirements

Python 3.10 or later is recommended.

Install the required Python packages with:

```bash
pip install -r requirements.txt
```

## Input format

A list of strings, each representing a clonotype, in the format:

CDR3a+Va;CDR3b+Vb

Up to 2 chains per locus can be included:

TRA1:TRA2;TRB1:TRB2

## Usage

```python
import tcr_epitope_3db

# Import databases
dbs = tcr_epitope_3db.import_3db(
    "path/to/vdjdb.tsv",
    "path/to/iedb.csv",
    "path/to/mcpas.csv"
)

# Match clonotypes
matches = tcr_epitope_3db.match_3db(
    clonotypes,
    *dbs,
    max_distance=0
)
```

`max_distance` controls the maximum Levenshtein distance allowed for CDR3 matching; it defaults to 0 for exact matching, while values greater than 0 allow approximate matches.
