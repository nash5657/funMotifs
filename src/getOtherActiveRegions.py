"""
Created on 05 Oct 2022

@author: Mark Melzer
"""

'''
Extract other active regions from bed files
'''

import os
from getInactivePromoters import get_cell_info_for_regions


def getBed(coordinates_input_file, input_file, output_file):
    oligos = open(coordinates_input_file, 'r').readlines()
    oligos_dict = {}
    for o in oligos:
        so = o.strip().split('\t')
        if so[3] not in oligos_dict.keys():
            oligos_dict[so[3]] = [so[0], so[1], so[2]]

    lines = open(input_file, 'r').readlines()
    out_file = open(output_file, 'w')

    for o in lines:
        so = o.strip().split(' ')
        if so[3] == 'ref' and (float(so[8]) >= 2 or float(so[13]) >= 2):
            out_file.write(
                '\t'.join(oligos_dict[so[1]]) + '\t' + so[0] + '\t' + str(max(float(so[5]), float(so[10]))) + '\n')

    return output_file


def get_other_active_regions(sys_args, db_name, col_names):
    Other_active_regions_file = sys_args[-4]
    if not os.path.exists(Other_active_regions_file):
        Other_active_regions_file = getBed(coordinates_input_file=sys_args[-6], input_file=sys_args[-5],
                                           output_file=Other_active_regions_file)
    cells_to_extract_info_from = sys_args[-3].split(',')
    regions_df_output_file = sys_args[-2]
    Other_active_regions_df = get_cell_info_for_regions(
        Other_active_regions_file, db_name=db_name, cells=cells_to_extract_info_from, assays=['all'],
        sep='\t', report_cols_from_file=True, cols_indices_to_report_from_file=[4],
        cols_names_to_report_from_file=['Activity_Score'], df_output_file=regions_df_output_file,
        region_name_index=3, region_strand_index=None, region_score_index=4, motif_score_index=4,
        max_number_motifs_to_report=1,
        min_dist_from_region_start=False, min_dist_from_region_center=True, max_motif_score=False,
        max_region_score=False,
        col_names=col_names)
    if not Other_active_regions_df.empty:
        Other_active_regions_df.to_csv(Other_active_regions_df + '.tsv', sep='\t')
    else:
        print("Is empty")


if __name__ == '__main__':
    db_name = "funmotifsdb"
    args = ['../datafiles/TrainingSets/Tewhey_Cell2016_20150102_Geuv.HepG2_SigRegs.bed', 'HepG2,Liver',
            '../datafiles/TrainingSets/HepG2Motifs_Tewhy.df', 'other_active_region']
    # prom_unactive_genes_start_index_params = 0
    col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
    get_other_active_regions(args, db_name, col_names)
    print("Done.")
    # TODO: write unit test
