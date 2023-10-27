"""
Created on Nov 28, 2022

@author: Mark Melzer

Receives the computed regression coefficients and a database with annotated motifs.
Output: List of Motifs that are functional according to the following conditions:
    - DHSs must be present
    - TF, belonging to motif, must be expressed
    - Functional if binding evidence from matching TFs
    - Non-functional when evidence of non-binding of matching TF
    - If no TF ChIP-seq data available: required significant functionality score
        --> significance based on distribution of scores from functional motifs (DHS and matching TF-peak) in
            myeloid tissue
"""
import psycopg2
import pandas as pd
import numpy as np
from DBUtilities import get_number_of_motifs, add_column_to_tissue_table, update_db_value, open_connection


# TODO: check why significance based on myeloid tissue and not tissue-wise?

def DHS_present(motif_info) -> bool:
    # TODO: check if dnase_seq is right variable and its possible values
    """
    Function that returns a boolean value to indicate whether DHSs is present for a motif
    """
    if motif_info == 'NaN' or motif_info == 'NO':
        return False
    elif motif_info == 'YES':
        return True
    else:
        print("DHS:", motif_info)
        raise Exception("Unexpected DHS format")


def TFs_expressed(motif_info) -> bool:
    """
    Function that returns a boolean value to indicate whether the TF belonging to motifs is expressed
    """
    if motif_info == 'NaN' or motif_info <= 0.0:
        return False

    return True


def TF_ChIP_seq_data_available(motif_info) -> bool:
    """
    Function that returns a boolean value to indicate whether TF ChIP-seq data is available for a motif
    """
    if motif_info == 'NaN':
        return False

    return True


def TF_binding_evidence(motif_info) -> bool:
    """
    Function that returns a boolean value to indicate whether there is evidence for binding of a matching TF to a motif
    """
    # TODO: find out what entry in table would be instead of checking all possibilities
    if motif_info in ['True', 'true', 'TRUE', 't', 'T', True, 1, 'YES', 'yes', 'Yes', 'y', 'Y']:
        return True

    return False


def compute_functionality_score(motif: object, params, weighted_variable: list, tissue: str, db_name: str,
                                db_user_name: str) -> float:
    """
    Function that computes the functionality score of a motif
    """
    try:
        score = params['intercept']
    except:
        score = 0
        pass

    for var in weighted_variable:
        # TODO: check type sensitivity
        if motif[var] == "NaN":
            # TODO: is there a better way to do this? --> values should be imputed before this step already
            continue
        # float removes information from data frame that is not necessary (like data type)
        score += float(motif[var]) * float(params[var])

    update_db_value(score, motif['mid'], 'funScore', tissue, db_name, db_user_name)
    return score


def get_significance_cutoff(funMotifs: list, nonFunMotifs: list, params, weighted_variable: list,
                            tissue, db_name, db_user_name) -> float:
    """
    Function that determines the cutoff value for a significant functionality score
    It is assumed that the funScore is normally distributed for functional and non-functional motifs.
    As cutoff, the interception between both distribution functions (density) is computed.
    """
    # get functional motifs in data frame
    funMotifs_df = get_motif_data_for_tissue(tissue, '*', db_name, db_user_name, mid=funMotifs, df=True)
    nonFunMotifs_df = get_motif_data_for_tissue(tissue, '*', db_name, db_user_name, mid=nonFunMotifs, df=True)

    funScores = []
    nonFunScores = []

    # TODO: number of motifs here must be sufficiently high

    # get the functionality scores for the functional motifs
    for idx, motif in funMotifs_df.iterrows():
        funScores.append(compute_functionality_score(motif, params, weighted_variable, tissue, db_name, db_user_name))

    # get the functionality scores for the not functional motifs
    for idx, motif in nonFunMotifs_df.iterrows():
        nonFunScores.append(compute_functionality_score(motif, params, weighted_variable, tissue, db_name,
                                                        db_user_name))

    # compute mean and std of scores
    mu2 = np.mean(funScores)
    sig2 = np.std(funScores)
    mu1 = np.mean(nonFunScores)
    sig1 = np.std(nonFunScores)

    # compute the significance cutoff (intersection of distributions of scores of functional vs. non-functional motifs)
    up = sig1 * sig2 * np.sqrt(
        2 * (sig2 ** 2 - sig1 ** 2) * np.log(sig2 / sig1) + (mu1 - mu2) ** 2) - sig1 ** 2 * mu2 + sig2 ** 2 * mu1
    low = sig2 ** 2 - sig1 ** 2

    return up / low


def make_string_from_list(lst):
    """This function makes a string out of a list of names, to use it in a sql statement"""
    sql = ""
    if type(lst) is str:
        return lst
    elif type(lst) is list:
        for l in lst:
            assert type(l) is str
            sql = sql + l + ", "
    sql = sql[:-2]

    return sql


def get_motif_data_for_tissue(tissue, columns, db_name, db_user_name, db_host_name, mid=None, df: bool = False) -> object:
    """
    Function that given a tissue returns the data of a database for this tissue as data frame
    """
    # TODO: or only return one motif at a time
    # establish connection
    #conn = psycopg2.connect(database=db_name, user=db_user_name)
    conn = open_connection(db_name, db_user_name, db_host_name)
    cur = conn.cursor()

    col = make_string_from_list(columns)

    # create sql command
    sql = f'''SELECT {col} FROM {tissue}'''

    if mid is not None:
        if type(mid) is int:
            sql = sql + f" WHERE mid={mid}"
        elif type(mid) is list:
            rows = "("
            for id in mid:
                assert type(id) is int
                rows = rows + str(id) + ", "
            sql = sql + f" WHERE mid IN {rows[:-2]})"
    # get data from data base
    cur.execute(sql)
    if df:
        result = pd.read_sql_query(sql, conn)
    else:
        result = cur.fetchall()

    # close connection
    conn.close()
    print("result of get motif data from tissue", tissue, "for motif ", mid, "  is: \n", result)
    return result


def get_functional_motifs(params, tissue, weighted_variables: list, db_name: str, db_user_name: str, db_host_name: str) -> list:
    """
    Function that returns functional motifs (mid value) based on annotated motifs and the regression weights
    """
    # add two columns to tissue tables to save functionality score and whether the motif is functional
    # TODO: check that the right values is added to the right columns
    add_column_to_tissue_table(tissue, db_name, db_user_name, db_host_name, col_name="funScore", col_type="real")
    add_column_to_tissue_table(tissue, db_name, db_user_name, db_host_name, col_name="functionality", col_type="text")

    # lists to store mids of motifs
    funMotifs = []
    nonFunMotifs = []
    compute_score = []

    # get number of motifs
    # TODO: remove hard-coded motifs
    num_motifs = get_number_of_motifs("motifs", db_name, db_user_name, db_host_name)

    # loop over motifs
    for motif in range(1, num_motifs + 1):
        if get_motif_data_for_tissue(tissue, 'dnase__seq', db_name, db_user_name, db_host_name, motif) == []:
            print("empty")
            continue
        if not DHS_present(get_motif_data_for_tissue(tissue, 'indexdhs', db_name, db_user_name, db_host_name, motif)[0][0]) or not \
                TFs_expressed(get_motif_data_for_tissue(tissue, 'tfexpr', db_name, db_user_name, db_host_name, motif)[0][0]):
            nonFunMotifs.append(motif)
        elif TF_ChIP_seq_data_available(get_motif_data_for_tissue(tissue, 'tfbinding', db_name, db_user_name, db_host_name, motif)
                                        [0][0]):
            if TF_binding_evidence(get_motif_data_for_tissue(tissue, 'tfbinding', db_name, db_user_name, db_host_name, motif)[0][0]):
                funMotifs.append(motif)
                update_db_value('YES', motif, 'functionality', tissue, db_name, db_user_name)
                update_db_value('YES', 'motifs', 'functionality', tissue, db_name, db_user_name)

            else:
                nonFunMotifs.append(motif)
                update_db_value('NO', motif, 'functionality', tissue, db_name, db_user_name)
        else:
            compute_score.append(motif)

    # TODO: if not done tissue-wise: move this to parent function and pass as argument
    # get significance cutoff for functionality score
    cutoff = get_significance_cutoff(funMotifs, nonFunMotifs, params, weighted_variables, tissue, db_name,
                                     db_user_name)
    for mid in compute_score:
        motif = get_motif_data_for_tissue(tissue, '*', db_name, db_user_name, db_host_name, mid)
        # TODO: save functionality score to database
        if compute_functionality_score(motif, params, weighted_variables, tissue, db_name, db_user_name) > cutoff:
            funMotifs.append(motif['mid'])
            update_db_value('YES', motif, 'functionality', tissue, db_name, db_user_name)
            update_db_value('YES', 'motifs', 'functionality', tissue, db_name, db_user_name)
        else:
            update_db_value('NO', motif, 'functionality', tissue, db_name, db_user_name)

    return funMotifs


def get_functional_motifs_per_tissue(params, tissues: list, weighted_variables=None,
                                     db_name: str = "funmotifsdb", db_user_name: str = "", db_host_name: str = "") -> dict:
    """
    Function that returns functional motifs per tissue based on annotated motifs and the regression weights
    """
    if weighted_variables is None:
        weighted_variables = ['ChromHMM'.lower(), 'DNase__seq'.lower(), 'FANTOM'.lower(), 'NumOtherTFBinding'.lower(),
                              'RepliDomain'.lower(), 'TFBinding'.lower(), 'TFExpr'.lower(), 'score'.lower(),
                              'footprints'.lower(), 'cCRE'.lower(), 'IndexDHS'.lower(), 'RegElem'.lower(), 'intercept']
    # assert that at least one tissue is given
    assert len(tissues) > 0

    add_column_to_tissue_table('motifs', db_name, db_user_name, db_host_name, col_name="functionality", col_type="text")

    # TODO: assertion that the column names of the tissue tables (except mid) corresponds to the weighted_variable list

    funMotifs_per_tissue = {}

    for tissue in tissues:
        print("tissue: ", tissue)
        print(funMotifs_per_tissue)
        funMotifs_per_tissue[tissue] = get_functional_motifs(params, tissue, weighted_variables, db_name, db_user_name, db_host_name)

    return funMotifs_per_tissue
