import json
import re
import argparse
import pickle

import pandas as pd
import numpy as np
import Levenshtein



def create_tcr_query_df(clonotypes):
    """
    Convert a list of clonotype strings into an expanded query dataframe.

    Each clonotype is formatted as TRA;TRB, with multiple chains within
    a locus separated by ':' and CDR3/V gene pairs separated by '+'.

    Multiple chains within a locus are expanded into all possible
    TRA x TRB combinations.

    Returns
    -------
    pd.DataFrame
        Query dataframe with clonotype and chain information.
    """

    rows = []

    for clonotype_id, clonotype in enumerate(clonotypes):

        # Split TRA and TRB
        loci = clonotype.split(';')

        tra = None
        trb = None

        for locus in loci:
            if not locus:
                continue

            chains = locus.split(':')

            # Infer locus from V gene
            first_v = chains[0].split('+')[1]

            if first_v.startswith('TRA'):
                tra = chains
            elif first_v.startswith('TRB'):
                trb = chains

        # Parse chains within each locus
        def parse_chains(chains):
            if chains is None:
                return [(pd.NA, pd.NA)]

            parsed = []

            for chain in chains:
                cdr3, v_gene = chain.split('+')
                parsed.append((cdr3, v_gene))

            return parsed

        tra_chains = parse_chains(tra)
        trb_chains = parse_chains(trb)

        # Generate all possible TRA x TRB combinations
        query_number = 0

        for tra_cdr3, tra_v in tra_chains:
            for trb_cdr3, trb_v in trb_chains:

                rows.append({
                    'clonotype_id': clonotype_id,
                    'query_id': f'{clonotype_id}.{query_number}',
                    'tra_cdr3': tra_cdr3,
                    'tra_v': tra_v,
                    'trb_cdr3': trb_cdr3,
                    'trb_v': trb_v
                })

                query_number += 1

    return pd.DataFrame(rows)


def import_iedb(path):

    iedb = pd.read_csv(
        path,
        skiprows=1,
        low_memory=False
    )

    def clean_text(x):

        if pd.isna(x):
            return ''

        return re.sub(r'\s+', '', str(x))

    def normalize_v_genes(curated, calculated):

        values = (
            clean_text(curated)
            + ','
            + clean_text(calculated)
        )

        return sorted(
            set(
                x
                for x in values.split(',')
                if x != ''
            )
        )

    def normalize_cdr3(curated, calculated):

        values = (
            clean_text(curated)
            + ','
            + clean_text(calculated)
        )

        return sorted(
            set(
                x if not x.startswith('C') else x[1:-1]
                for x in values.split(',')
                if x != ''
            )
        )

    iedb['V Gene'] = iedb.apply(
        lambda row: normalize_v_genes(
            row['Curated V Gene'],
            row['Calculated V Gene']
        ),
        axis=1
    )

    iedb['CDR3'] = iedb.apply(
        lambda row: normalize_cdr3(
            row['CDR3 Curated'],
            row['CDR3 Calculated']
        ),
        axis=1
    )

    iedb['V Gene.1'] = iedb.apply(
        lambda row: normalize_v_genes(
            row['Curated V Gene.1'],
            row['Calculated V Gene.1']
        ),
        axis=1
    )

    iedb['CDR3.1'] = iedb.apply(
        lambda row: normalize_cdr3(
            row['CDR3 Curated.1'],
            row['CDR3 Calculated.1']
        ),
        axis=1
    )

    chain_a = pd.DataFrame({
        'iedb_row': iedb.index,
        'cdr3': iedb['CDR3'],
        'v_gene': iedb['V Gene'],
        'locus': 'TRA'
    })

    chain_b = pd.DataFrame({
        'iedb_row': iedb.index,
        'cdr3': iedb['CDR3.1'],
        'v_gene': iedb['V Gene.1'],
        'locus': 'TRB'
    })

    chain_a = chain_a[
        chain_a['cdr3'].apply(len) > 0
    ].copy()

    chain_b = chain_b[
        chain_b['cdr3'].apply(len) > 0
    ].copy()

    chain_a = chain_a.explode(
        'cdr3',
        ignore_index=True
    )

    chain_b = chain_b.explode(
        'cdr3',
        ignore_index=True
    )

    iedb_chains = pd.concat(
        [
            chain_a,
            chain_b
        ],
        ignore_index=True
    )

    return iedb, iedb_chains

def match_iedb(query_df, iedb, iedb_chains, max_distance=0):

    query = query_df.copy()
    query_columns = query_df.columns.tolist()

    # Normalize query CDR3s
    def normalize_cdr3(x):

        if pd.isna(x):
            return pd.NA

        x = str(x).replace(' ', '')

        if x.startswith('C'):
            return x[1:-1]

        return x

    query['_tra_cdr3'] = query['tra_cdr3'].apply(
        normalize_cdr3
    )

    query['_trb_cdr3'] = query['trb_cdr3'].apply(
        normalize_cdr3
    )

    # IEDB information to retain
    iedb_info_columns = [
        'IEDB Receptor ID',
        'Name',
        'Source Molecule',
        'Source Organism',
        'MHC Allele Names'
    ]

    iedb_info = iedb[
        iedb_info_columns
    ].copy()

    iedb_info.index.name = 'iedb_id'
    iedb_info = iedb_info.reset_index()

    # EXACT MATCHES
    exact_matches = []

    for locus, query_cdr3_col, query_v_col in [
        ('TRA', '_tra_cdr3', 'tra_v'),
        ('TRB', '_trb_cdr3', 'trb_v')
    ]:

        q = query[
            query[query_cdr3_col].notna()
        ][[
            'query_id',
            'clonotype_id',
            query_cdr3_col,
            query_v_col
        ]].copy()

        q.rename(
            columns={
                query_cdr3_col: 'query_cdr3',
                query_v_col: 'query_v'
            },
            inplace=True
        )

        db = iedb_chains[
            iedb_chains['locus'] == locus
        ][[
            'iedb_row',
            'cdr3',
            'v_gene'
        ]].copy()

        db = db[
            db['cdr3'].notna()
        ]

        matches = q.merge(
            db,
            left_on='query_cdr3',
            right_on='cdr3',
            how='inner'
        )

        if matches.empty:
            continue

        matches['match_locus'] = locus
        matches['match_distance'] = 0

        matches.rename(
            columns={
                'iedb_row': 'iedb_id'
            },
            inplace=True
        )

        if locus == 'TRA':

            matches['match_cdr3a'] = matches['cdr3']
            matches['match_va'] = matches['v_gene']

            matches['match_cdr3b'] = pd.NA
            matches['match_vb'] = pd.NA

        else:

            matches['match_cdr3a'] = pd.NA
            matches['match_va'] = pd.NA

            matches['match_cdr3b'] = matches['cdr3']
            matches['match_vb'] = matches['v_gene']

        exact_matches.append(
            matches[
                [
                    'query_id',
                    'clonotype_id',
                    'match_locus',
                    'match_distance',
                    'iedb_id',
                    'match_cdr3a',
                    'match_va',
                    'match_cdr3b',
                    'match_vb'
                ]
            ]
        )

    # DISTANCE MATCHES
    distance_matches = []

    if max_distance > 0:

        for locus, query_cdr3_col, query_v_col, db_cdr3_col in [
            ('TRA', '_tra_cdr3', 'tra_v', 'CDR3'),
            ('TRB', '_trb_cdr3', 'trb_v', 'CDR3.1')
        ]:

            q = query[
                query[query_cdr3_col].notna()
            ][[
                'query_id',
                'clonotype_id',
                query_cdr3_col,
                query_v_col
            ]].copy()

            q.rename(
                columns={
                    query_cdr3_col: 'query_cdr3',
                    query_v_col: 'query_v'
                },
                inplace=True
            )

            # IEDB stores CDR3s as lists
            cdr3_db = sorted(
                set(
                    x
                    for y in iedb[db_cdr3_col].tolist()
                    for x in y
                )
            )

            cdr3_query = sorted(
                q['query_cdr3'].unique().tolist()
            )

            if not cdr3_query or not cdr3_db:
                continue

            # Levenshtein distance matrix
            dist = np.empty(
                (
                    len(cdr3_query),
                    len(cdr3_db)
                ),
                dtype=np.int16
            )

            for i, query_cdr3 in enumerate(cdr3_query):

                for j, db_cdr3 in enumerate(cdr3_db):

                    dist[i, j] = Levenshtein.distance(
                        query_cdr3,
                        db_cdr3
                    )

            # Keep non-exact matches within max_distance
            rows, cols = np.where(
                (dist > 0) &
                (dist <= max_distance)
            )

            if len(rows) == 0:
                continue

            distance_pairs = pd.DataFrame({
                'query_cdr3': [
                    cdr3_query[i]
                    for i in rows
                ],
                'cdr3': [
                    cdr3_db[j]
                    for j in cols
                ],
                'match_distance': dist[rows, cols]
            })

            matches = q.merge(
                distance_pairs,
                on='query_cdr3',
                how='inner'
            )

            db = iedb_chains[
                iedb_chains['locus'] == locus
            ][[
                'iedb_row',
                'cdr3',
                'v_gene'
            ]].copy()

            db = db[
                db['cdr3'].notna()
            ]

            matches = matches.merge(
                db,
                on='cdr3',
                how='inner'
            )

            matches.rename(
                columns={
                    'iedb_row': 'iedb_id'
                },
                inplace=True
            )

            matches['match_locus'] = locus

            if locus == 'TRA':

                matches['match_cdr3a'] = matches['cdr3']
                matches['match_va'] = matches['v_gene']

                matches['match_cdr3b'] = pd.NA
                matches['match_vb'] = pd.NA

            else:

                matches['match_cdr3a'] = pd.NA
                matches['match_va'] = pd.NA

                matches['match_cdr3b'] = matches['cdr3']
                matches['match_vb'] = matches['v_gene']

            distance_matches.append(
                matches[
                    [
                        'query_id',
                        'clonotype_id',
                        'match_locus',
                        'match_distance',
                        'iedb_id',
                        'match_cdr3a',
                        'match_va',
                        'match_cdr3b',
                        'match_vb'
                    ]
                ]
            )

    # Combine exact + distance matches
    pieces = exact_matches + distance_matches

    if not pieces:

        return pd.DataFrame(
            columns=query_columns + [
                'match_locus',
                'match_distance',
                'match_paired_distance',
                'iedb_id',
                'match_cdr3a',
                'match_va',
                'match_cdr3b',
                'match_vb'
            ] + iedb_info_columns
        )

    matches = pd.concat(
        pieces,
        ignore_index=True
    )

    # Add query information
    matches = matches.merge(
        query[query_columns],
        on=[
            'query_id',
            'clonotype_id'
        ],
        how='left'
    )

    # Calculate paired-chain distance
    def calculate_paired_distance(row):

        if row['match_locus'] == 'TRA':

            query_paired = normalize_cdr3(
                row['trb_cdr3']
            )

            db_paired = row['match_cdr3b']

        else:

            query_paired = normalize_cdr3(
                row['tra_cdr3']
            )

            db_paired = row['match_cdr3a']

        if pd.isna(query_paired) or pd.isna(db_paired):
            return pd.NA

        return Levenshtein.distance(
            query_paired,
            db_paired
        )

    matches['match_paired_distance'] = matches.apply(
        calculate_paired_distance,
        axis=1
    )

    # Add selected IEDB information
    matches = matches.merge(
        iedb_info,
        on='iedb_id',
        how='left'
    )

    # Final column order
    matches = matches[
        query_columns
        + [
            'match_locus',
            'match_distance',
            'match_paired_distance',
            'iedb_id',
            'match_cdr3a',
            'match_va',
            'match_cdr3b',
            'match_vb'
        ]
        + iedb_info_columns
    ]

    return matches


def import_mcpas(path):

    mcpas = pd.read_csv(path,
                        low_memory=False)

    def remove_spaces(x):

        if pd.isna(x):
            return np.nan

        return re.sub(r'\s+', '', str(x))

    def normalize_cdr3(x):

        x = remove_spaces(x)

        if pd.isna(x) or x == '':
            return np.nan

        # Remove external C residues when present
        if x.startswith('C'):
            return x[1:-1]

        return x

    # Normalize CDR3 alpha
    mcpas['CDR3.alpha.aa.inner'] = (
        mcpas['CDR3.alpha.aa']
        .apply(normalize_cdr3)
    )

    # Normalize CDR3 beta
    mcpas['CDR3.beta.aa.inner'] = (
        mcpas['CDR3.beta.aa']
        .apply(normalize_cdr3)
    )

    # Clean TRAV
    #
    # Original McPAS annotation is preserved, including alleles.
    mcpas['TRAV'] = (
        mcpas['TRAV']
        .apply(remove_spaces)
    )

    # Clean TRBV
    #
    # Original McPAS annotation is preserved, including alleles.
    mcpas['TRBV'] = (
        mcpas['TRBV']
        .apply(remove_spaces)
    )

    # Create chain A table
    chain_a = pd.DataFrame({
        'mcpas_row': mcpas.index,
        'cdr3': mcpas['CDR3.alpha.aa.inner'],
        'v_gene': mcpas['TRAV'],
        'locus': 'TRA'
    })

    # Create chain B table
    chain_b = pd.DataFrame({
        'mcpas_row': mcpas.index,
        'cdr3': mcpas['CDR3.beta.aa.inner'],
        'v_gene': mcpas['TRBV'],
        'locus': 'TRB'
    })

    # Remove rows without CDR3
    chain_a = chain_a[
        chain_a['cdr3'].notna()
    ].copy()

    chain_b = chain_b[
        chain_b['cdr3'].notna()
    ].copy()

    # Explode CDR3 values
    chain_a = chain_a.explode(
        'cdr3',
        ignore_index=True
    )

    chain_b = chain_b.explode(
        'cdr3',
        ignore_index=True
    )

    # Combine
    mcpas_chains = pd.concat(
        [
            chain_a,
            chain_b
        ],
        ignore_index=True
    )

    return mcpas, mcpas_chains


def match_mcpas(query_df, mcpas, mcpas_chains, max_distance=0):

    query = query_df.copy()
    query_columns = query_df.columns.tolist()

    def normalize_cdr3(x):

        if pd.isna(x):
            return pd.NA

        x = str(x).replace(' ', '')

        if x.startswith('C'):
            return x[1:-1]

        return x

    query['_tra_cdr3'] = query['tra_cdr3'].apply(
        normalize_cdr3
    )

    query['_trb_cdr3'] = query['trb_cdr3'].apply(
        normalize_cdr3
    )

    mcpas_info_columns = [
        'Species',
        'Category',
        'Pathology',
        'Antigen.protein',
        'Epitope.peptide',
        'Epitope.ID',
        'MHC',
        'Tissue',
        'T.Cell.Type',
        'PubMed.ID',
        'Remarks'
    ]

    mcpas_info = mcpas[
        mcpas_info_columns
    ].copy()

    mcpas_info.index.name = 'mcpas_id'
    mcpas_info = mcpas_info.reset_index()

    exact_matches = []

    for locus, query_cdr3_col, query_v_col in [
        ('TRA', '_tra_cdr3', 'tra_v'),
        ('TRB', '_trb_cdr3', 'trb_v')
    ]:

        q = query[
            query[query_cdr3_col].notna()
        ][[
            'query_id',
            'clonotype_id',
            query_cdr3_col,
            query_v_col
        ]].copy()

        q.rename(
            columns={
                query_cdr3_col: 'query_cdr3',
                query_v_col: 'query_v'
            },
            inplace=True
        )

        db = mcpas_chains[
            mcpas_chains['locus'] == locus
        ][[
            'mcpas_row',
            'cdr3',
            'v_gene'
        ]].copy()

        db = db[
            db['cdr3'].notna()
        ]

        matches = q.merge(
            db,
            left_on='query_cdr3',
            right_on='cdr3',
            how='inner'
        )

        if matches.empty:
            continue

        matches['match_locus'] = locus
        matches['match_distance'] = 0

        matches.rename(
            columns={
                'mcpas_row': 'mcpas_id'
            },
            inplace=True
        )

        if locus == 'TRA':

            matches['match_cdr3a'] = matches['cdr3']
            matches['match_va'] = matches['v_gene']

            matches['match_cdr3b'] = pd.NA
            matches['match_vb'] = pd.NA

        else:

            matches['match_cdr3a'] = pd.NA
            matches['match_va'] = pd.NA

            matches['match_cdr3b'] = matches['cdr3']
            matches['match_vb'] = matches['v_gene']

        exact_matches.append(
            matches[
                [
                    'query_id',
                    'clonotype_id',
                    'match_locus',
                    'match_distance',
                    'mcpas_id',
                    'match_cdr3a',
                    'match_va',
                    'match_cdr3b',
                    'match_vb'
                ]
            ]
        )

    distance_matches = []

    if max_distance > 0:

        for locus, query_cdr3_col, query_v_col in [
            ('TRA', '_tra_cdr3', 'tra_v'),
            ('TRB', '_trb_cdr3', 'trb_v')
        ]:

            q = query[
                query[query_cdr3_col].notna()
            ][[
                'query_id',
                'clonotype_id',
                query_cdr3_col,
                query_v_col
            ]].copy()

            q.rename(
                columns={
                    query_cdr3_col: 'query_cdr3',
                    query_v_col: 'query_v'
                },
                inplace=True
            )

            cdr3_query = sorted(
                q['query_cdr3'].unique().tolist()
            )

            cdr3_db = sorted(
                mcpas_chains.loc[
                    mcpas_chains['locus'] == locus,
                    'cdr3'
                ].dropna().unique().tolist()
            )

            if not cdr3_query or not cdr3_db:
                continue

            dist = np.empty(
                (
                    len(cdr3_query),
                    len(cdr3_db)
                ),
                dtype=np.int16
            )

            for i, query_cdr3 in enumerate(cdr3_query):

                for j, db_cdr3 in enumerate(cdr3_db):

                    dist[i, j] = Levenshtein.distance(
                        query_cdr3,
                        db_cdr3
                    )

            rows, cols = np.where(
                (dist > 0) &
                (dist <= max_distance)
            )

            if len(rows) == 0:
                continue

            distance_pairs = pd.DataFrame({
                'query_cdr3': [
                    cdr3_query[i]
                    for i in rows
                ],
                'cdr3': [
                    cdr3_db[j]
                    for j in cols
                ],
                'match_distance': dist[rows, cols]
            })

            matches = q.merge(
                distance_pairs,
                on='query_cdr3',
                how='inner'
            )

            db = mcpas_chains[
                mcpas_chains['locus'] == locus
            ][[
                'mcpas_row',
                'cdr3',
                'v_gene'
            ]].copy()

            db = db[
                db['cdr3'].notna()
            ]

            matches = matches.merge(
                db,
                on='cdr3',
                how='inner'
            )

            matches.rename(
                columns={
                    'mcpas_row': 'mcpas_id'
                },
                inplace=True
            )

            matches['match_locus'] = locus

            if locus == 'TRA':

                matches['match_cdr3a'] = matches['cdr3']
                matches['match_va'] = matches['v_gene']

                matches['match_cdr3b'] = pd.NA
                matches['match_vb'] = pd.NA

            else:

                matches['match_cdr3a'] = pd.NA
                matches['match_va'] = pd.NA

                matches['match_cdr3b'] = matches['cdr3']
                matches['match_vb'] = matches['v_gene']

            distance_matches.append(
                matches[
                    [
                        'query_id',
                        'clonotype_id',
                        'match_locus',
                        'match_distance',
                        'mcpas_id',
                        'match_cdr3a',
                        'match_va',
                        'match_cdr3b',
                        'match_vb'
                    ]
                ]
            )

    pieces = exact_matches + distance_matches

    if not pieces:

        return pd.DataFrame(
            columns=query_columns + [
                'match_locus',
                'match_distance',
                'match_paired_distance',
                'mcpas_id',
                'match_cdr3a',
                'match_va',
                'match_cdr3b',
                'match_vb'
            ] + mcpas_info_columns
        )

    matches = pd.concat(
        pieces,
        ignore_index=True
    )

    matches = matches.merge(
        query[query_columns],
        on=[
            'query_id',
            'clonotype_id'
        ],
        how='left'
    )

    def calculate_paired_distance(row):

        if row['match_locus'] == 'TRA':

            query_paired = normalize_cdr3(
                row['trb_cdr3']
            )

            db_paired = row['match_cdr3b']

        else:

            query_paired = normalize_cdr3(
                row['tra_cdr3']
            )

            db_paired = row['match_cdr3a']

        if pd.isna(query_paired) or pd.isna(db_paired):
            return pd.NA

        return Levenshtein.distance(
            query_paired,
            db_paired
        )

    matches['match_paired_distance'] = matches.apply(
        calculate_paired_distance,
        axis=1
    )

    matches = matches.merge(
        mcpas_info,
        on='mcpas_id',
        how='left'
    )

    matches = matches[
        query_columns
        + [
            'match_locus',
            'match_distance',
            'match_paired_distance',
            'mcpas_id',
            'match_cdr3a',
            'match_va',
            'match_cdr3b',
            'match_vb'
        ]
        + mcpas_info_columns
    ]

    return matches


def import_vdjdb(path):

    vdjdb = pd.read_csv(
        path,
        sep='\t',
        low_memory=False
    )

    def remove_spaces(x):

        if pd.isna(x):
            return np.nan

        return re.sub(r'\s+', '', str(x))

    def normalize_cdr3(x):

        if pd.isna(x):
            return np.nan

        x = remove_spaces(x)

        if x == '':
            return np.nan

        if x.startswith('C'):
            return x[1:-1]

        return x

    # Clean all string columns
    #
    # This removes spaces, including non-breaking spaces and
    # other whitespace characters, while leaving the content
    # otherwise unchanged.
    for col in vdjdb.columns:

        if vdjdb[col].dtype == 'object':
            vdjdb[col] = vdjdb[col].apply(remove_spaces)

    # Extract CDR3fix information
    #
    # CDR3fix is stored as a JSON string containing:
    #   cdr3
    #   cdr3_old
    #
    # We retain both representations.
    def get_cdr3fix_values(x):

        if pd.isna(x):
            return []

        try:
            fix = json.loads(x)
        except (json.JSONDecodeError, TypeError):
            return []

        values = []

        for key in ['cdr3', 'cdr3_old']:

            value = fix.get(key)

            if pd.notna(value):
                value = normalize_cdr3(value)

                if pd.notna(value):
                    values.append(value)

        return sorted(set(values))

    vdjdb['CDR3.inner'] = (
        vdjdb['CDR3fix']
        .apply(get_cdr3fix_values)
    )

    # Create VDJdb chain table
    #
    # VDJdb has one chain per row.
    #
    # Gene determines the locus:
    #   TRB â†’ TRB
    #   TRA â†’ TRA
    #   TRG â†’ TRG
    #   TRD â†’ TRD
    vdjdb_chains = pd.DataFrame({
        'vdjdb_row': vdjdb.index,
        'complex_id': vdjdb['complex.id'],
        'cdr3': vdjdb['CDR3.inner'],
        'v_gene': vdjdb['V'],
        'locus': vdjdb['Gene']
    })

    # Remove rows without CDR3 annotations
    vdjdb_chains = vdjdb_chains[
        vdjdb_chains['cdr3'].apply(len) > 0
    ].copy()

    # Explode CDR3 lists
    #
    # If CDR3 and CDR3_old differ, both become separate rows.
    # If they are identical, only one row is created.
    vdjdb_chains = vdjdb_chains.explode(
        'cdr3',
        ignore_index=True
    )

    return vdjdb, vdjdb_chains


def match_vdjdb(query_df, vdjdb, vdjdb_chains, max_distance=0):

    query = query_df.copy()
    query_columns = query_df.columns.tolist()

    def normalize_cdr3(x):

        if pd.isna(x):
            return pd.NA

        x = str(x).replace(' ', '')

        if x.startswith('C'):
            return x[1:-1]

        return x

    query['_tra_cdr3'] = query['tra_cdr3'].apply(
        normalize_cdr3
    )

    query['_trb_cdr3'] = query['trb_cdr3'].apply(
        normalize_cdr3
    )

    vdjdb_info_columns = [
        'Epitope',
        'Epitope gene',
        'Epitope species',
        'Species',
        'MHC A',
        'MHC B',
        'MHC class',
        'Reference'
    ]

    vdjdb_info = vdjdb[
        vdjdb_info_columns
    ].copy()

    vdjdb_info.index.name = 'vdjdb_id'

    vdjdb_info = vdjdb_info.reset_index()

    exact_matches = []

    for locus, query_cdr3_col in [
        ('TRA', '_tra_cdr3'),
        ('TRB', '_trb_cdr3')
    ]:

        q = query[
            query[query_cdr3_col].notna()
        ][[
            'query_id',
            'clonotype_id',
            query_cdr3_col
        ]].copy()

        q.rename(
            columns={
                query_cdr3_col: 'query_cdr3'
            },
            inplace=True
        )

        db = vdjdb_chains[
            vdjdb_chains['locus'] == locus
        ][[
            'vdjdb_row',
            'complex_id',
            'cdr3',
            'v_gene'
        ]].copy()

        db = db[
            db['cdr3'].notna()
        ]

        matches = q.merge(
            db,
            left_on='query_cdr3',
            right_on='cdr3',
            how='inner'
        )

        if matches.empty:
            continue

        matches['match_locus'] = locus
        matches['match_distance'] = 0

        matches.rename(
            columns={
                'vdjdb_row': 'vdjdb_id'
            },
            inplace=True
        )

        if locus == 'TRA':

            matches['match_cdr3a'] = matches['cdr3']
            matches['match_va'] = matches['v_gene']

            matches['match_cdr3b'] = pd.NA
            matches['match_vb'] = pd.NA

        else:

            matches['match_cdr3a'] = pd.NA
            matches['match_va'] = pd.NA

            matches['match_cdr3b'] = matches['cdr3']
            matches['match_vb'] = matches['v_gene']

        exact_matches.append(
            matches[
                [
                    'query_id',
                    'clonotype_id',
                    'match_locus',
                    'match_distance',
                    'vdjdb_id',
                    'match_cdr3a',
                    'match_va',
                    'match_cdr3b',
                    'match_vb'
                ]
            ]
        )

    distance_matches = []

    if max_distance > 0:

        for locus, query_cdr3_col in [
            ('TRA', '_tra_cdr3'),
            ('TRB', '_trb_cdr3')
        ]:

            q = query[
                query[query_cdr3_col].notna()
            ][[
                'query_id',
                'clonotype_id',
                query_cdr3_col
            ]].copy()

            q.rename(
                columns={
                    query_cdr3_col: 'query_cdr3'
                },
                inplace=True
            )

            db = vdjdb_chains[
                (vdjdb_chains['locus'] == locus)
                &
                (vdjdb_chains['cdr3'].notna())
            ][[
                'vdjdb_row',
                'complex_id',
                'cdr3',
                'v_gene'
            ]].copy()

            cdr3_db = sorted(
                db['cdr3'].unique().tolist()
            )

            cdr3_query = sorted(
                q['query_cdr3'].unique().tolist()
            )

            if not cdr3_query or not cdr3_db:
                continue

            dist = np.empty(
                (
                    len(cdr3_query),
                    len(cdr3_db)
                ),
                dtype=np.int16
            )

            for i, query_cdr3 in enumerate(cdr3_query):

                for j, db_cdr3 in enumerate(cdr3_db):

                    dist[i, j] = Levenshtein.distance(
                        query_cdr3,
                        db_cdr3
                    )

            rows, cols = np.where(
                (dist > 0)
                &
                (dist <= max_distance)
            )

            if len(rows) == 0:
                continue

            distance_pairs = pd.DataFrame({
                'query_cdr3': [
                    cdr3_query[i]
                    for i in rows
                ],
                'cdr3': [
                    cdr3_db[j]
                    for j in cols
                ],
                'match_distance': dist[rows, cols]
            })

            matches = q.merge(
                distance_pairs,
                on='query_cdr3',
                how='inner'
            )

            matches = matches.merge(
                db,
                on='cdr3',
                how='inner'
            )

            matches['match_locus'] = locus

            matches.rename(
                columns={
                    'vdjdb_row': 'vdjdb_id'
                },
                inplace=True
            )

            if locus == 'TRA':

                matches['match_cdr3a'] = matches['cdr3']
                matches['match_va'] = matches['v_gene']

                matches['match_cdr3b'] = pd.NA
                matches['match_vb'] = pd.NA

            else:

                matches['match_cdr3a'] = pd.NA
                matches['match_va'] = pd.NA

                matches['match_cdr3b'] = matches['cdr3']
                matches['match_vb'] = matches['v_gene']

            distance_matches.append(
                matches[
                    [
                        'query_id',
                        'clonotype_id',
                        'match_locus',
                        'match_distance',
                        'vdjdb_id',
                        'match_cdr3a',
                        'match_va',
                        'match_cdr3b',
                        'match_vb'
                    ]
                ]
            )

    pieces = exact_matches + distance_matches

    if not pieces:

        return pd.DataFrame(
            columns=query_columns + [
                'match_locus',
                'match_distance',
                'match_paired_distance',
                'vdjdb_id',
                'match_cdr3a',
                'match_va',
                'match_cdr3b',
                'match_vb'
            ] + vdjdb_info_columns
        )

    matches = pd.concat(
        pieces,
        ignore_index=True
    )

    # Add VDJdb complex ID
    #
    # This is done after exact + distance matches have been
    # combined.
    matches = matches.merge(
        vdjdb_chains[
            [
                'vdjdb_row',
                'complex_id'
            ]
        ].drop_duplicates(),
        left_on='vdjdb_id',
        right_on='vdjdb_row',
        how='left'
    )

    matches.drop(
        columns='vdjdb_row',
        inplace=True
    )

    # Annotate paired chain
    #
    # Only nonzero complex IDs represent paired chains.
    paired_db = vdjdb_chains[
        vdjdb_chains['complex_id'] != 0
    ][[
        'complex_id',
        'locus',
        'cdr3',
        'v_gene'
    ]].copy()

    # TRA chains
    tra_pairs = paired_db[
        paired_db['locus'] == 'TRA'
    ][[
        'complex_id',
        'cdr3',
        'v_gene'
    ]].copy()

    tra_pairs.rename(
        columns={
            'cdr3': 'paired_cdr3a',
            'v_gene': 'paired_va'
        },
        inplace=True
    )

    # TRB chains
    trb_pairs = paired_db[
        paired_db['locus'] == 'TRB'
    ][[
        'complex_id',
        'cdr3',
        'v_gene'
    ]].copy()

    trb_pairs.rename(
        columns={
            'cdr3': 'paired_cdr3b',
            'v_gene': 'paired_vb'
        },
        inplace=True
    )

    # Merge paired chains
    matches = matches.merge(
        tra_pairs,
        on='complex_id',
        how='left'
    )

    matches = matches.merge(
        trb_pairs,
        on='complex_id',
        how='left'
    )

    # Fill paired chain information
    tra_mask = matches['match_locus'] == 'TRA'
    trb_mask = matches['match_locus'] == 'TRB'

    matches.loc[
        tra_mask,
        'match_cdr3b'
    ] = matches.loc[
        tra_mask,
        'paired_cdr3b'
    ]

    matches.loc[
        tra_mask,
        'match_vb'
    ] = matches.loc[
        tra_mask,
        'paired_vb'
    ]

    matches.loc[
        trb_mask,
        'match_cdr3a'
    ] = matches.loc[
        trb_mask,
        'paired_cdr3a'
    ]

    matches.loc[
        trb_mask,
        'match_va'
    ] = matches.loc[
        trb_mask,
        'paired_va'
    ]

    matches.drop(
        columns=[
            'paired_cdr3a',
            'paired_va',
            'paired_cdr3b',
            'paired_vb'
        ],
        inplace=True
    )

    # Add query information
    #
    # Needed before calculating paired-chain distance.
    matches = matches.merge(
        query[query_columns],
        on=[
            'query_id',
            'clonotype_id'
        ],
        how='left'
    )

    # Calculate paired-chain distance
    #
    # The primary match distance is between the matched chain
    # and its corresponding query chain.
    #
    # This distance is between the query's OTHER chain and the
    # paired chain found in VDJdb.
    def calculate_paired_distance(row):

        if row['match_locus'] == 'TRA':

            query_paired = normalize_cdr3(
                row['trb_cdr3']
            )

            db_paired = row['match_cdr3b']

        else:

            query_paired = normalize_cdr3(
                row['tra_cdr3']
            )

            db_paired = row['match_cdr3a']

        if pd.isna(query_paired) or pd.isna(db_paired):
            return pd.NA

        return Levenshtein.distance(
            query_paired,
            db_paired
        )

    matches['match_paired_distance'] = matches.apply(
        calculate_paired_distance,
        axis=1
    )

    # Add selected VDJdb information
    matches = matches.merge(
        vdjdb_info,
        on='vdjdb_id',
        how='left'
    )

    # Final column order
    matches = matches[
        query_columns
        + [
            'match_locus',
            'match_distance',
            'match_paired_distance',
            'vdjdb_id',
            'match_cdr3a',
            'match_va',
            'match_cdr3b',
            'match_vb'
        ]
        + vdjdb_info_columns
    ]

    return matches


def annotate_mhc(x):

    if pd.isna(x):
        return pd.Series(
            [np.nan, np.nan],
            index=['mhc_class', 'species']
        )

    s = str(x).upper().replace('\xa0', ' ').strip()

    parts = [
        p.strip()
        for p in re.split(r'[|/]', s)
        if p.strip()
    ]

    classes = []
    species = []

    for p in parts:

        # Human MHC-I
        if (
            re.search(r'\bHLA[- ]?[A-C E-G](?:\*|[0-9]|CW|\b)', p)
            or re.search(
                r'\b(?:HUMAN )?CD1[ABCDE]?\b|\b(?:HUMAN )?MR1\b',
                p
            )
        ):
            classes.append('MHCI')
            species.append('Human')

        # Human MHC-II
        elif re.search(r'\bHLA[- ]?(?:DR|DQ|DP)', p):
            classes.append('MHCII')
            species.append('Human')

        # Mouse MHC-I
        elif re.search(r'\bH-?2[- ]?[KDLQ][A-Z0-9]*\b', p):
            classes.append('MHCI')
            species.append('Mouse')

        # Mouse MHC-II
        elif re.search(
            r'\bI-[AE][A-Z0-9]*\b|\bH-?2[- ]?I[AE][A-Z0-9]*\b',
            p
        ):
            classes.append('MHCII')
            species.append('Mouse')

        # Generic class labels
        elif re.search(r'\bHLA CLASS I\b', p):
            classes.append('MHCI')
            species.append('Human')

        elif re.search(r'\bHLA CLASS II\b', p):
            classes.append('MHCII')
            species.append('Human')

        elif re.search(r'\bH-?2 CLASS I\b', p):
            classes.append('MHCI')
            species.append('Mouse')

        elif re.search(r'\bH-?2 CLASS II\b', p):
            classes.append('MHCII')
            species.append('Mouse')

    mhc_class = (
        classes[0]
        if classes and all(c == classes[0] for c in classes)
        else np.nan
    )

    species_out = (
        species[0]
        if species and all(s == species[0] for s in species)
        else np.nan
    )

    return pd.Series(
        [mhc_class, species_out],
        index=['mhc_class', 'species']
    )


def combine_matches(
    matches_iedb,
    matches_mcpas,
    matches_vdjdb
):

    common_columns = [
        'clonotype_id',
        'query_id',
        'tra_cdr3',
        'tra_v',
        'trb_cdr3',
        'trb_v',
        'match_locus',
        'match_distance',
        'match_paired_distance',
        'match_cdr3a',
        'match_va',
        'match_cdr3b',
        'match_vb'
    ]

    # IEDB
    iedb = matches_iedb[common_columns].copy()

    iedb['match_db'] = 'iedb'
    iedb['db_row'] = matches_iedb['iedb_id']
    iedb['epitope'] = matches_iedb['Name']
    iedb['antigen'] = matches_iedb['Source Molecule']
    iedb['source'] = matches_iedb['Source Organism']
    iedb['mhc'] = matches_iedb['MHC Allele Names']

    if not iedb.empty:
        iedb[['mhc_class', 'species']] = (
            iedb['mhc'].apply(annotate_mhc)
        )
    else:
        iedb['mhc_class'] = pd.Series(dtype='object')
        iedb['species'] = pd.Series(dtype='object')

    # McPAS
    mcpas = matches_mcpas[common_columns].copy()

    mcpas['match_db'] = 'mcpas'
    mcpas['db_row'] = matches_mcpas['mcpas_id']
    mcpas['epitope'] = matches_mcpas['Epitope.peptide']
    mcpas['antigen'] = matches_mcpas['Antigen.protein']
    mcpas['source'] = matches_mcpas['Pathology']
    mcpas['mhc'] = matches_mcpas['MHC']

    if not mcpas.empty:
        mcpas['mhc_class'] = (
            mcpas['mhc']
            .apply(lambda x: annotate_mhc(x)['mhc_class'])
        )
    else:
        mcpas['mhc_class'] = pd.Series(dtype='object')
    mcpas['species'] = matches_mcpas['Species']

    # McPAS already has species
    mcpas['species'] = matches_mcpas['Species']

    # VDJdb
    vdjdb = matches_vdjdb[common_columns].copy()

    vdjdb['match_db'] = 'vdjdb'
    vdjdb['db_row'] = matches_vdjdb['vdjdb_id']
    vdjdb['epitope'] = matches_vdjdb['Epitope']
    vdjdb['antigen'] = matches_vdjdb['Epitope gene']
    vdjdb['source'] = matches_vdjdb['Epitope species']

    # Combine MHC A and MHC B
    vdjdb['mhc'] = (
        matches_vdjdb[['MHC A', 'MHC B']]
        .apply(
            lambda row: '|'.join(
                str(x)
                for x in row
                if pd.notna(x) and str(x).strip()
            ),
            axis=1
        )
        .replace('', np.nan)
    )

    # VDJdb already provides these
    vdjdb['mhc_class'] = matches_vdjdb['MHC class']
    vdjdb['species'] = matches_vdjdb['Species']

    columns = common_columns + [
        'match_db',
        'db_row',
        'epitope',
        'antigen',
        'source',
        'mhc',
        'mhc_class',
        'species'
    ]

    return pd.concat(
        [
            iedb[columns],
            mcpas[columns],
            vdjdb[columns]
        ],
        ignore_index=True
    ).sort_values(by='query_id')


def import_3db(path_vdjdb, path_iedb, path_mcpas):
    """
    Import and preprocess the VDJdb, IEDB, and McPAS-TCR databases.

    Parameters
    ----------
    path_vdjdb : str
        Path to the VDJdb file.
    path_iedb : str
        Path to the IEDB file.
    path_mcpas : str
        Path to the McPAS-TCR file.

    Returns
    -------
    tuple
        Imported and preprocessed database objects and chain tables.
    """

    vdjdb, vdjdb_chains = import_vdjdb(path_vdjdb)
    iedb, iedb_chains = import_iedb(path_iedb)
    mcpas, mcpas_chains = import_mcpas(path_mcpas)

    return (
        vdjdb,
        vdjdb_chains,
        iedb,
        iedb_chains,
        mcpas,
        mcpas_chains
    )


def match_3db(
    clonotypes,
    vdjdb, vdjdb_chains,
    iedb, iedb_chains,
    mcpas, mcpas_chains,
    max_distance=0
):
    """
    Match TCR clonotypes against VDJdb, IEDB, and McPAS-TCR.

    Parameters
    ----------
    clonotypes : list of str
        TCR clonotypes in the expected input format.
    vdjdb, vdjdb_chains : pandas.DataFrame
        Imported VDJdb data.
    iedb, iedb_chains : pandas.DataFrame
        Imported IEDB data.
    mcpas, mcpas_chains : pandas.DataFrame
        Imported McPAS-TCR data.
    max_distance : int, default=0
        Maximum Levenshtein distance allowed for CDR3 matching.

    Returns
    -------
    pandas.DataFrame
        Table containing TCR-epitope matches.
    """

    query_df = create_tcr_query_df(clonotypes)

    matches_vdjdb = match_vdjdb(
        query_df,
        vdjdb,
        vdjdb_chains,
        max_distance=max_distance
    )

    matches_iedb = match_iedb(
        query_df,
        iedb,
        iedb_chains,
        max_distance=max_distance
    )

    matches_mcpas = match_mcpas(
        query_df,
        mcpas,
        mcpas_chains,
        max_distance=max_distance
    )

    matches_comb = combine_matches(
        matches_vdjdb=matches_vdjdb,
        matches_iedb=matches_iedb,
        matches_mcpas=matches_mcpas
    )

    return matches_comb


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Match TCR clonotypes against VDJdb, IEDB, and McPAS-TCR.'
    )

    subparsers = parser.add_subparsers(
        dest='command',
        required=True
    )

    # Import command
    import_parser = subparsers.add_parser(
        'import',
        help='Import and preprocess the databases and save them as a pickle.'
    )

    import_parser.add_argument(
        '--vdjdb',
        required=True,
        help='Path to the VDJdb database file.'
    )

    import_parser.add_argument(
        '--iedb',
        required=True,
        help='Path to the IEDB database file.'
    )

    import_parser.add_argument(
        '--mcpas',
        required=True,
        help='Path to the McPAS-TCR database file.'
    )

    import_parser.add_argument(
        '-o',
        '--output',
        default='tcr_databases.pkl',
        help='Output pickle file (default: tcr_databases.pkl).'
    )

    # Match command
    match_parser = subparsers.add_parser(
        'match',
        help='Match clonotypes against previously imported databases.'
    )

    match_parser.add_argument(
        'clonotypes',
        help='Path to a text file containing one clonotype per line.'
    )

    match_parser.add_argument(
        '--db',
        required=True,
        help='Path to the pickle file containing imported databases.'
    )

    match_parser.add_argument(
        '-d',
        '--distance',
        type=int,
        default=0,
        help='Maximum Levenshtein distance for CDR3 matching (default: 0).'
    )

    match_parser.add_argument(
        '-o',
        '--output',
        default='tcr_matches.csv',
        help='Output CSV file (default: tcr_matches.csv).'
    )

    args = parser.parse_args()

    # Import databases
    if args.command == 'import':

        dbs = import_3db(
            args.vdjdb,
            args.iedb,
            args.mcpas
        )

        with open(args.output, 'wb') as f:
            pickle.dump(
                dbs,
                f
            )

        print(
            f'Database import complete. '
            f'Processed databases saved to {args.output}'
        )

    # Match clonotypes
    elif args.command == 'match':

        with open(args.clonotypes) as f:
            clonotypes = [
                line.strip()
                for line in f
                if line.strip()
            ]

        with open(args.db, 'rb') as f:
            dbs = pickle.load(f)

        matches = match_3db(
            clonotypes,
            *dbs,
            max_distance=args.distance
        )

        if matches.empty:
            print(
                'No matches found in VDJdb, IEDB, or McPAS-TCR. '
                'The output is empty.'
            )

        matches.to_csv(
            args.output,
            index=False
        )

        print(
            f'Matching complete. Results saved to {args.output}'
        )
