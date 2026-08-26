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

## Output

## Output

The function returns a pandas DataFrame containing the identified TCR–epitope matches across VDJdb, IEDB, and McPAS-TCR, including query clonotype information, matched CDR3 and V-gene sequences, matching distances, epitope and antigen information, and MHC annotations.

## Databases
The tool currently supports matching against VDJdb, IEDB, and McPAS-TCR. Database files must be downloaded separately and provided to the import function. Each database is processed using database-specific import and matching logic to account for differences in their data structures.

## Citation
If you use this tool in your research, please cite the original databases used.
