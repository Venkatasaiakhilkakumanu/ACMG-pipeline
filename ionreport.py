#!/usr/bin/env python3
# ionreport.py

from urllib.request import urlopen
import fileinput
import sys
from io import StringIO
import csv
import pandas as pd
import numpy as np
import argparse
import re
import io


def main():
    report()


def ion():
   df = pd.read_csv("/media/bioserve/CGI9/6-6-23/PB3135_ICP_germline-full.tsv", sep='\t', header=0, skiprows=2)
   df.rename(columns={df.columns[0]: 'locus'}, inplace=True)
   cs = df.iloc[:, [0, 1, 2, 4, 9, 13, 18, 19, 24, 25, 31, 41]]
   cs[['chr', 'pos']] = cs['locus'].str.split(':', expand=True)
   cs['chr'] = cs['chr'].str.replace('chr', '')
   cs[['genotype', 'alt']] = cs['genotype'].str.split('/', expand=True)

   cs['BA1'] = ''
   cs['BS1'] = ''
   cs['BP4'] = ''
   cs['BP4'] = ''
   cs['BP6'] = ''
   cs['PVS1'] = ''
   cs['PP3'] = ''
   cs['PP5'] = ''
   cs['PM2'] = ''

   cs['Classification'] = ''

   file = "/media/bioserve/CGI9/Akhil/pb3137.txt"
   cs.to_csv(file, sep='\t', index=False)
   return(cs)


def report():
    cs = ion()

    ba = pd.read_csv("/media/bioserve/CGI9/Akhil/pb3137.txt", sep='\t', header=0)
    ba = pd.DataFrame(ba)
    print(ba)

    def get_lowest_value(row):
        if isinstance(row, str) and ':' in row:
            values = [float(value) for value in row.split(':')]
            return min(values)
        else:
            return float(row)


    ba['maf'] = ba['maf'].apply(get_lowest_value)
    ba['BA1'] = ba['maf'].apply(lambda x: 1 if x > 0.05 else 0)

    #BA1
    #ba['BA1'] = ba['maf'].apply(lambda x: 1 if x > 0.05 else 0)


    #PVS1
    patho = ["frameshiftDeletion", "frameshiftInsertion", "stoploss", "nnsense", "nonsense"]
    patho2 = "splice"
    ba['PVS1'] = 0
    ba.loc[ba['function'].isin(patho), 'PVS1'] = 1
    ba.loc[ba['location'].str.contains(patho2, flags=re.IGNORECASE, regex=True, na=False), 'PVS1'] = 1


    #PM2
    column_mapping = {
        '0': 'Combined_column',
        '1_x': 'gnomad',
        '1_y': 'allele_count',
        '2': 'dom_rec'
    }

    db = pd.read_csv('/media/bioserve/CGI9/Akhil/pm2_drc.txt', sep='\t', names=column_mapping.values(), header=0)
    #db_genome = pd.read_csv('/media/bioserve/CGI9/Akhil/Alleles/gnomad_genome_maf_combined.txt', sep='\t', names=['Combined_column', 'gnomad_genome'], header=0)
    ba['Combined_column'] = ba.iloc[:, [12, 13, 2, 14]].apply(lambda x: '_'.join(map(str, x)), axis=1)
    merged_df = pd.merge(ba, db, on='Combined_column', how='left')
    #merged_df_gene = pd.merge(ba, db_genome, on='Combined_column', how='left')
    ba['gnomad'] = merged_df['gnomad']
    #ba['gnomad_genome'] = merged_df_gene['gnomad_genome']
    ba['allele count'] = merged_df['allele_count']
    ba['dom_rec'] = merged_df['dom_rec']
    ba['gnomad'] = ba['gnomad'].fillna(np.nan)  # Replace empty cells with NaN
    ba['gnomad'] = pd.to_numeric(ba['gnomad'], errors='coerce')  # Convert column to numeric values
    ba['PM2'] = np.where(pd.isnull(ba['gnomad']) & (ba['coverage'] > 20), 1, 0)
    # #ba['PM2'] = ba['gnomad'].apply(lambda x: 0 if pd.notnull(x) and ba['coverage'] < 78.1 else 0)

    #PP3
    ba['phylop'] = ba['phylop'].astype(str)
    ba['PP3'] = 0
    # Iterate over each row in the DataFrame
    for index, row in ba.iterrows():
        phylop_value = row['phylop']
        # Skip rows with multiple values separated by ","
        if ',' in phylop_value:
            continue
        # Convert the phylop value to float
        try:
            phylop_value = float(phylop_value)
        except ValueError:
            continue
        # Assign 1 if phylop value is greater than 0, otherwise assign 0
        if phylop_value > 7.52:
            ba.at[index, 'PP3'] = 1

    #PP5
    clinvar = ["Likely Pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic"]
    ba["PP5"] = 0
    ba.loc[ba['clinvar'].isin(clinvar), 'PP5'] = 1

   #BS1
    for index, row in ba.iterrows():
        maf_value = row['maf']
        ba1_value = row['BA1']
        # Check if BA1 value is 1, if so, skip assigning a value to BS1
        if ba1_value == 1:
            continue
        # Check if maf value is greater than 0.005
        if maf_value > 0.005:
            ba.at[index, 'BS1'] = 1

    #BP4
    ba['phylop'] = ba['phylop'].astype(str)
    ba['BP4'] = 0
    # Iterate over each row in the DataFrame
    for index, row in ba.iterrows():
        phylop_value = row['phylop']
        # Skip rows with multiple values separated by ","
        if ',' in phylop_value:
            continue
        # Convert the phylop value to float
        try:
            phylop_value = float(phylop_value)
        except ValueError:
            continue
        # Assign 1 if phylop value is greater than 0, otherwise assign 0
        if phylop_value < 3.58:
            ba.at[index, 'BP4'] = 1

    #BP6
    clinvar_b = ["Likely Benign", "Benign"]
    ba['BP6'] = 0
    ba.loc[ba['clinvar'].isin(clinvar_b), 'BP6'] = 1

    #Classification
    ba['Classification'] = ba['BA1'].apply(lambda x: 'Benign' if x == 1 else '')
    like = (ba['BS1'] == 1) & (ba['BP4'] == 1)
    ba.loc[like, 'Classification'] = 'Likely Benign'
    patho5 = (ba['PVS1'] == 1) & (ba['PM2'] == 1) & (ba['PP5'] == 1)
    patho3 = (ba['PVS1'] == 1) & (ba['PM2'] == 1) & (ba['PP3'] == 1)
    likelypatho = (ba['PVS1'] == 1) & (ba['PM2'] == 1) & (ba['PP5'] == 0)
    ba.loc[likelypatho, 'Classification'] = 'Likely Pathogenic'
    ba.loc[patho5, 'Classification'] = 'Pathogenic'
    ba.loc[patho3, 'Classification'] = 'Pathogenic'

    file1 = "/media/bioserve/CGI9/Akhil/pb3137_classified.txt"
    ba.to_csv(file1, sep='\t', index=False)


if __name__ == "__main__":
    main()
