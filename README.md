# tcr_epitope_annotation

Utilities for matching T-cell receptor (TCR) clonotypes against three
TCR–epitope databases:

- VDJdb
- IEDB
- McPAS-TCR

The code supports exact CDR3 matching and optional Levenshtein distance
matching.

## Requirements

Python 3.10 or later is recommended.

Install the required Python packages with:

```bash
pip install -r requirements.txt
```

## Input format

First, databases have to be imported by passing their local paths to the import command. Then, TCRs are passed to the match command as a list of strings, each representing a clonotype in the bioidentity format: CDR3+V.

Locus is inferred from the V gene name. Paired chains can be represented as: CDR3a+Va;CDR3b+Vb.

Up to 2 chains per locus can be included: TRA1:TRA2;TRB1:TRB2.

## Usage

1. Import databases
```bash
python tcr_epitope_3db.py import --vdjdb path/to/vdjdb.tsv --iedb path/to/iedb.csv --mcpas path/to/mcpas.csv
```

2. Match clonotypes
```bash
python tcr_epitope_3db.py match clonotypes.txt --db tcr_databases.pkl
```

`max_distance` controls the maximum Levenshtein distance allowed for CDR3 matching; it defaults to 0 for exact matching, while values greater than 0 allow approximate matches.

## Output

The function returns tcr_matches.csv containing the identified TCR-epitope matches across VDJdb, IEDB, and McPAS-TCR, including query clonotype information, matched CDR3 and V gene sequences, matching distances, epitope and antigen information, and MHC annotations.

## Databases
The tool currently supports matching against VDJdb, IEDB, and McPAS-TCR. Database files must be downloaded separately and provided to the import function. Each database is processed using database-specific import and matching logic to account for differences in their data structures. The databases can be downloaded from their respective websites:

- VDJdb: https://vdjdb.com/search (export as TSV with no filters set)
- IEDB: https://www.iedb.org/downloader.php?file_name=doc/tcr_full_v3.zip
- McPAS-TCR: https://friedmanlab.weizmann.ac.il/McPAS-TCR/

In this way, VDJdb should be TSV while IEDB and McPAS-TCR should be CSV.
## Citation
If you use this tool in your research, please cite this repository and the original databases used.
