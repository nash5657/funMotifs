"""
Created on Nov 29, 2022

@author: Mark Melzer

Input: list of functional motifs, database in which the annotated motifs are saved, file containing non-coding mutations
Output: subset of these mutations that occur in the functional motifs
"""

import pybedtools
import psycopg2
import pandas as pd
import os

"""
Idea: 

- get functional motifs in bed format
- get variants in bed format
- use bed tool intersect

BedFormat:
chr start_pos    end_pos  name (motif: motif name, variant: tumor-sample??)    strand   mid (for motifs)
"""


def get_Bedtool_from_dataframe(df: object):
    """
    Function that takes a dataframe and outputs a BedTool object
    """
    return pybedtools.BedTool.from_dataframe(df)


def get_BedTool_for_functional_motifs(funMotifs: dict, tissue: str, db_user_name: str, db_name: str):
    """
    Function that extracts the information for motifs functional in a specified tissue from a psql database
    Output: BedTool object created from a data frame
    """
    # establish connection
    conn = psycopg2.connect(
        database=db_name, user=db_user_name)

    # create list of motif ids that will be extracted from table
    motifs = "("
    for funMotif in funMotifs[tissue]:
        motifs = motifs + f"{funMotif}, "
    motifs = motifs[:-2] + ")"

    # create sql command
    sql = f"SELECT {'chr', 'motifstart', 'motifend', 'name', 'strand', 'mid'} FROM {'motifs'} WHERE mid IN {motifs}"

    # get data into data frame
    df = pd.read_sql_query(sql, conn)
    df.set_index('mid', inplace=True)

    # close connection
    conn.close()

    return get_Bedtool_from_dataframe(df)


def get_BedTool_from_variant_file(variant_file: str):
    """
    Function that extracts the essential information of a variant file and returns it as BedTool object
    """
    # TODO: create checks for file format
    os.system("awk '{print $2, $3, $4, $9}' " + variant_file + f" >{variant_file}_tmp")
    return pybedtools.BedTool(variant_file + "_tmp")


def overlap_variants_and_motifs(motifs, variants, output_file: str):
    """
    Function that overlaps variants and functional motifs
    """
    return motifs.intersect(variants, wo=True).saveas(output_file)


def find_funMotif_variants_in_tissue(funMotifs: dict, tissue: str, variant_BedTool: object, db_name: str,
                                     db_user_name: str, output_file: str):
    """
    Function that returns the overlaps between functional motifs of a tissue and variants
    """
    motif_BedTool = get_BedTool_for_functional_motifs(funMotifs, tissue, db_user_name, db_name)
    funMotif_variants = overlap_variants_and_motifs(motif_BedTool, variant_BedTool, output_file)
    return funMotif_variants


def find_funMotif_variants(funMotifs: dict, tissues: list, variant_file: str, db_user_name: str, output_file: str,
                           db_name: str = "funmotifsdb", force_overwrite: bool = False):
    """
    Function that returns the overlaps between functional motifs of all specified tissue and variants
    """
    variant_BedTool = get_BedTool_from_variant_file(variant_file)

    assert len(tissues) > 0
    for tissue in tissues:
        assert type(tissue) == str
        output_file_tissue = output_file + tissue
        if not os.path.exists(output_file_tissue) or force_overwrite:
            find_funMotif_variants_in_tissue(funMotifs, tissue, variant_BedTool, db_name, db_user_name,
                                             output_file_tissue)

    # TODO: check where return values are needed (save_as, intersect, etc.)
    return
