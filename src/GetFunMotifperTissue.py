"""
Created on Nov 28, 2022

@author: Mark Melzer

Receives the computed regression coefficients and a database with annotated motifs.
Output: List of Motifs that are functional according to the regression.
Condition:
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

# TODO: check why significance based on myeloid tissue and not tissue-wise?

def DHS_present(motif: object) -> bool:
    # TODO: check if dnase_seq is right variable
    """
    Function that returns a boolean value to indicate whether DHSs is present for a motif
    """
    if motif['dnase__seq'] == 'NaN' or motif['dnase_seq'] <= 0.0:
        return False

    return True


def TFs_expressed(motif: object) -> bool:
    """
    Function that returns a boolean value to indicate whether the TF belonging to motifs is expressed
    """
    if motif['tfexpr'] == 'NaN' or motif['tfexpr'] <= 0.0:
        return False

    return True


def TF_ChIP_seq_data_available(motif: object) -> bool:
    # TODO: check if tfbinding is right variable
    """
    Function that returns a boolean value to indicate whether TF ChIP-seq data is available for a motif
    """
    if motif['tfbinding'] == 'NaN':
        return False

    return True


def TF_binding_evidence(motif: object) -> bool:
    # TODO: check if tfbinding is right variable
    """
    Function that returns a boolean value to indicate whether there is evidence for binding of a matching TF to a motif
    """
    # TODO: find out what entrance in table would be instead of checking all possiblilities
    if motif['tfbinding'] in ['True', 'true', 'TRUE', 't', 'T', True, 1, 'YES', 'yes', 'Yes', 'y', 'Y']:
        return True

    return False


def compute_functionality_score(motif: object, params, weighted_variable: list) -> float:
    """
    Function that computes the functionality score of a motif
    """
    # TODO: check what x0 is called
    score = params[0]
    for var in weighted_variable:
        score += motif[var] * params[var]

    return score


def get_significance_cutoff() -> float:
    """
    Function that determines the cutoff value for a significant functionality score
    """
    return


def get_motif_data_for_tissue(tissue, db_name, db_user_name) -> object:
    """
    Function that given a tissue returns the data of a database for this tissue as data frame
    """
    # TODO: or only return one motif at a time
    # establish connection
    conn = psycopg2.connect(
        database=db_name, user=db_user_name)

    # create sql command
    sql = f'''SELECT * from {tissue}'''

    # get data into data frame
    df = pd.read_sql_query(sql, conn)
    df.reset_index()

    # close connection
    conn.close()

    return df


def get_functional_motifs(params, tissue, weighted_variables: list, db_name: str, db_user_name: str) -> list:
    """
    Function that returns functional motifs (mid value) based on annotated motifs and the regression weights
    """
    # get motifs from database (tissue tables)
    # TODO: might be more efficient to extract column by column (#motifs should be same in each tissue) and compute
    motifs = get_motif_data_for_tissue()
    funMotifs = []
    compute_score = []

    # loop over rows of the annotated motif data frame
    for idx, motif in motifs.iterrows():
        # TODO: check condition priority after next meeting
        if DHS_present(motif):
            funMotifs.append(motif['mid'])
        elif not TFs_expressed(motif):
            continue
        elif TF_ChIP_seq_data_available():
            if TF_binding_evidence():
                funMotifs.append(motif['mid'])
        else:
            compute_score.append(motif)

    # TODO: if not done tissue-wise: move this to parent function and pass as argument
    # get significance cutoff for functionality score
    cutoff = get_significance_cutoff(funMotifs)
    for motif in compute_score:
        if compute_functionality_score(motif) >= cutoff:
            funMotifs.append(motif['mid'])

    return funMotifs


def get_functional_motifs_per_tissue(params, tissues: list, weighted_variables: list = ['all'],
                                     db_name: str = "funmotifsdb", db_user_name: str = "") -> list:
    """
    Function that returns functional motifs per tissue based on annotated motifs and the regression weights
    """
    if weighted_variables == ['all']:
        # noinspection PyStatementEffect
        weighted_variables == ['ChromHMM'.lower(), 'DNase__seq'.lower(), 'FANTOM'.lower(), 'NumOtherTFBinding'.lower(),
                               'RepliDomain'.lower(), 'TFBinding'.lower(), 'TFExpr'.lower(), 'score'.lower(),
                               'footprints'.lower(), 'cCRE'.lower(), 'IndexDHS'.lower(), 'RegElem'.lower()]
    # assert that at least one tissue is given
    assert len(tissues) > 0

    # TODO: assertion that the column names of the tissue tables (except mid) corresponds to the weighted_variable list

    funMotifs_per_tissue = [{}] * len(tissues)

    for index, tissue in enumerate(tissues):
        funMotifs_per_tissue[index][tissue] = get_functional_motifs(params, tissue, weighted_variables, db_name,
                                                                    db_user_name)

    return funMotifs_per_tissue
