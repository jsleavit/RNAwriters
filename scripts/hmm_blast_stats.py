#!/usr/bin/env python

import sys
import os
import glob 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os,sys


# gtRNAdb info
genomeInfo = pd.read_csv('/projects/lowelab/users/jsleavit/git_repos/data/Gtrnadb_prokaryotic_info_with_full_taxonomy_updated_asm_accessions.tsv', sep = '\t', index_col = 0)

# blast results
blastResults = pd.read_csv('/projects/lowelab/users/jsleavit/git_repos/data/uniprot_modomics_mod_enzyme_hits-dedup.tsv', sep = '\t')

# build dict of the curated modomics database
modomicsModDb = pd.read_csv('/projects/lowelab/users/jsleavit/git_repos/data/modomics_cog_combined_database.tsv', sep = '\t')
uniprotModDict = {'data':[]}
for i in range(len(modomicsModDb)):
    uniprot = modomicsModDb['uniprot'].values[i]
    cog = modomicsModDb['Unnamed: 0'].values[i]
    symbol = modomicsModDb['symbol'].values[i]
    if type(modomicsModDb['modifications'].values[i]) != float:
        try: 
             # if there are multiple uniprot ids for a single COG id write then all out as their own entry in the dictionary
            uniprot = uniprot.split(';')
            for u in uniprot:
                uniprotModDict['data'].append({'COG':cog,'uniprot':u, 'symbol':symbol, 'modifications':list(set([x.split(':')[-1] for x in modomicsModDb['modifications'].values[i].split(';') if 't' in x]))})
        except:
            uniprotModDict['data'].append({'COG':cog,'uniprot':uniprot, 'symbol':symbol, 'modifications':list(set([x.split(':')[-1] for x in modomicsModDb['modifications'].values[i].split(';') if 't' in x]))})


# get the predicted mod enzymes and read in the gzipped file
predictedModEnzymes = pd.read_csv('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/test_hmm_blast_results/hmm_blast_results.gz', sep = '\t', compression = 'gzip')

''' split the query and clean up the column annotations and formatting '''

# split the query column into multiple columns by the '|' delimiter
split_cols = predictedModEnzymes.iloc[:, 0].str.split('|', expand=True)
max_cols = split_cols.shape[1] + predictedModEnzymes.shape[1]
print(max_cols, predictedModEnzymes.shape[1])
if max_cols > predictedModEnzymes.shape[1]:
    for i in range(max_cols):

        # add empty columns to the dataframe
        predictedModEnzymes[f'col_{predictedModEnzymes.shape[1]+i}'] = ''

# add the original column to the dataframes
predictedModEnzymes.iloc[:, split_cols.shape[1]:max_cols] = predictedModEnzymes.iloc[:, 1:predictedModEnzymes.shape[1]]

# add the split columns to the dataframe
predictedModEnzymes.iloc[:, 0:split_cols.shape[1]] = split_cols
predModEnzSplit = predictedModEnzymes 

# reasign the column names
predModEnzSplit.rename(columns = {'query': 'gene_acc','target':'hmm_cog','percent_identity':'hmm_evalue',
                                   '3': 'gb_acc', '4': 'phylum', '5': 'query_product', '7': 'target', '8': 'percent_identity',
                                   '9':'length', 'e_value': 'mismatch', 'bit_score':'gapopen', 'target_product':'qstart', 'col_13': 'qend',
                                   'col_15': 'sstart', 'col_17': 'send', 'col_19': 'blast_evalue', 'col_21': 'bitscore', 'col_23': 'target_product'}, inplace = True)


# if column 6 is not empty, concat it to column 5, then drop column 6 (some products have a '|' in them), and the empty columns
predModEnzSplit['query_product'] = predModEnzSplit['query_product'].fillna('') + predModEnzSplit['6'].fillna('')
predModEnzSplit.drop(columns = ['6'], inplace = True)
predModEnzSplit.drop(columns = [x for x in predModEnzSplit.columns if 'col' in x], inplace = True)


# add columns from the uniprot curated modomics database (this is pretty slow, maybe try implementing vectorization)
columns_to_add = {'hmm_symbol':'symbol', 'blast_symbol':'symbol', 'modifications':'modifications','blast_cog':'COG','sum_eval':'sum_eval'}
for c in columns_to_add:

    if c not in predModEnzSplit.columns:

        if 'sum_eval' in c:
            predModEnzSplit[c] = predModEnzSplit['hmm_evalue'].astype(float) + predModEnzSplit['blast_evalue'].astype(float)

        if 'blast' in c:

            access = columns_to_add[c]
            predModEnzSplit[c] = predModEnzSplit['target'].apply(lambda x: [y[access] for y in uniprotModDict['data'] if y['uniprot'] == x])
            predModEnzSplit[c] = predModEnzSplit[c].apply(lambda x: x[0] if len(x) > 0 else 'NA')

        if 'hmm' in c:
            access = columns_to_add[c]
            predModEnzSplit[c] = predModEnzSplit['hmm_cog'].apply(lambda x: [y[access] for y in uniprotModDict['data'] if y['COG'] == x])
            predModEnzSplit[c] = predModEnzSplit[c].apply(lambda x: x[0] if len(x) > 0 else 'NA')

        if 'modifications' in c:
            access = columns_to_add[c]
            predModEnzSplit[c] = predModEnzSplit['target'].apply(lambda x: [y[access] for y in uniprotModDict['data'] if y['uniprot'] == x]) 
            predModEnzSplit[c] = predModEnzSplit[c].apply(lambda x: [n for y in x for n in y] if len(x) > 0 else 'NA')

# write out to a gz file for later use
predModEnzSplit.to_csv('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/hmm_blast_results_annotated.gz', sep = '\t', compression = 'gzip', index = False)

# subset for when hmm_symbol is equal to the blast_symbol
predModEnzSplit_sub = predModEnzSplit[predModEnzSplit['hmm_symbol'] == predModEnzSplit['blast_symbol']]

