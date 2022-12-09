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


def get_Bedtool_from_dataframe(df: object, output_file: str):
    """
    Function that takes a dataframe and outputs a BedTool object
    """
    pybedtools.BedTool.from_dataframe(df).saveas(
        output_file
    )
    return


def get_BedTool_for_functional_motifs(funMotifs: dict, tissue: str, db_user_name: str, db_name: str, output_file: str):
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

    # TODO: remove necessity of motifs, but somehow contain tissue-dependency
    # create sql command
    sql = f"SELECT chr, motifstart, motifend, name, strand, mid FROM motifs WHERE mid IN {motifs} " \
          f"AND functionality = 'YES'"

    # get data into data frame
    df = pd.read_sql_query(sql, conn)

    # close connection
    conn.close()

    get_Bedtool_from_dataframe(df, output_file)

    return


def get_BedTool_from_variant_file(variant_file: str):
    """
    Function that extracts the essential information of a variant file and returns it as BedTool object
    """
    # TODO: create checks for file format
    os.system(f"grep -v '#' {variant_file} " + "| awk -F '\t' '{print substr($2, 4), $3, $4, $9}' " +
              f" >{variant_file}_tmp")
    pybedtools.BedTool(variant_file + "_tmp").saveas(variant_file + ".bed_tmp")
    os.system(f"perl -p -i -e 's/ /\t/g' {variant_file}.bed_tmp")
    os.system(f"sed '1d' {variant_file}.bed_tmp > {variant_file}.bed")
    os.remove(f"{variant_file}.bed_tmp")
    return


def overlap_variants_and_motifs(motifs, variants, output_file: str):
    """
    Function that overlaps variants and functional motifs
    """
    # TODO: if necessary to make BedTool again, change architecture, figure out why one file okay and the other not
    mot = pybedtools.BedTool(motifs)
    mot.intersect(variants, wo=True, header=True).saveas(output_file)
    return


def find_funMotif_variants_in_tissue(funMotifs: dict, tissue: str, variant_BedTool_file: str, db_name: str,
                                     db_user_name: str, output_file: str, motif_BedTool_file: str):
    """
    Function that returns the overlaps between functional motifs of a tissue and variants
    """
    get_BedTool_for_functional_motifs(funMotifs, tissue, db_user_name, db_name, motif_BedTool_file)
    overlap_variants_and_motifs(motif_BedTool_file, variant_BedTool_file, output_file)
    return


def find_funMotif_variants(funMotifs: dict, tissues: list, variant_file: str, db_user_name: str, output_file: str,
                           db_name: str = "funmotifsdb", force_overwrite: bool = False):
    """
    Function that returns the overlaps between functional motifs of all specified tissue and variants
    """
    get_BedTool_from_variant_file(variant_file)
    motif_BedTool_file_short = "funMotif_in_"

    assert len(tissues) > 0
    for tissue in tissues:
        assert type(tissue) == str
        output_file_tissue = output_file + '_' + tissue
        if not os.path.exists(output_file_tissue) or force_overwrite:
            motif_BedTool_file = motif_BedTool_file_short + tissue + ".bed"
            find_funMotif_variants_in_tissue(funMotifs, tissue, variant_file + ".bed", db_name, db_user_name,
                                             output_file_tissue, motif_BedTool_file)

    # TODO: check where return values are needed (save_as, intersect, etc.)
    # TODO: can make it once for all functional motifs and prevent computing motifs multiple times
    return
