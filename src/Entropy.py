"""
Created on Aug 15, 2023

@author: Mark Melzer

Input: File containing mutation in functional motifs
Output: Change of entropy caused by this mutation
"""

import csv

import psycopg2


def mutation_entropy(row, db_name, db_user_name):
    """
    Computes the entropy for a mutation
    """
    # compute position of mutation in motif
    position_in_motif = compute_position_in_motif(row)

    # connect to database
    conn = psycopg2.connect(database=db_name, user=db_user_name)
    curs = conn.cursor()

    # create sql command
    query = f"SELECT freq FROM motifs_pfm WHERE position = {position_in_motif} AND name = {row['Name']} AND" \
                    f" allele = {row['Reference_Allele']}" \
                    f" - SELECT freq FROM motifs_pfm WHERE position = {position_in_motif} AND name = {row['Name']}" \
                    f" AND allele = {row['Mutated_Allele']}"
    curs.execute(query)

    return curs.fetchone()


def compute_position_in_motif(row):
    """
    Computes the position of the mutation in the motif
    """
    if row['Strand'] == "+":
        return row['Motif_Start'] - row['Variant_Start']
    else:
        return row['Variant_Start'] - row['Motif_End']


def compute_entropy(infile, outfile, db_name, db_user_name):
    """
    Computes entropy change caused by mutations
    """
    # read file
    with open(infile, 'r') as infile, open(outfile, 'w') as outfile:
        header = ["Chr", "Motif_Start", "Motif_End", "Name", "Strand", "mid", "Chr", "Variant_Start", "Variant_End",
                  "Variant_Classification", "Variant_Type", "Reference_Allele", "Mutated_Allele",
                  "Tumour_Sample_Barcode", "#Overlapping_Base_Pairs"]
        reader = csv.DictReader(infile, delimiter='\t', fieldnames=header)
        fieldnames = reader.fieldnames + ['Entropy']

        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            # set entropy to one for structural mutations
            if row['Variant_Type'] != "SNP":
                row['Entropy'] = 1
                writer.writerow(row)
            # compute entropy for SNVs
            else:
                row['Entropy'] = mutation_entropy(row, db_name, db_user_name)
                writer.writerow(row)

    return
