'''
Created on Dec 9, 2016

@author: husensofteng
'''
import sys, os
from pybedtools import BedTool
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import numpy as np

from GetDBData_TrainingSets import get_cell_info_for_motifs, get_cell_info_for_regions

def score_per_pos(MPRA_tiles_input_file, output_file, experiment_cell_name, tile_length = 145, region_pos_col = 0, region_chr_index = 3, cell_name_index = 0, region_center_index=4, region_name_sep='_', values_start_index = 1, sep = '\t'):
    
    with open(MPRA_tiles_input_file, 'r') as infile, open(output_file, 'w') as outfile:
        line = infile.readline()
        while line:
            out_lines = []
            sline = line.split()
            region_chr = sline[region_pos_col].split(region_name_sep)[region_chr_index]
            region_center = sline[region_pos_col].split(region_name_sep)[region_center_index]
            values = sline[values_start_index::]
            curr_pos = tile_length+2
            before_center = True
            for v in values:
                if curr_pos==0:
                    before_center = False
                if before_center:#set the start and end for each region
                    out_lines.append([region_chr, str(int(region_center)-(tile_length-curr_pos)), str(int(region_center)-(tile_length-curr_pos)), experiment_cell_name, v, str(curr_pos), sline[region_pos_col]])
                    curr_pos-=1
                else:
                    out_lines.append([region_chr, str(int(region_center)+curr_pos), str(int(region_center)+curr_pos), experiment_cell_name, v, str(curr_pos), sline[region_pos_col]])
                    curr_pos+=1
            line = infile.readline()
            for out_line in out_lines:
                outfile.write(sep.join(out_line) + '\n')
    return

def get_motif_scores(scores_per_bp_input_file, motifs_dir, motifs_scored_output_file, max_motif_score_per_region=False):
    
    regions = {}
    scores_per_bp_input_file_obj = BedTool(scores_per_bp_input_file)
    for motif_file in os.listdir(motifs_dir):
        motif_file_obj = BedTool(motifs_dir+'/'+motif_file)
        results = motif_file_obj.intersect(scores_per_bp_input_file_obj, wo=True).sort().groupby(g=[1,2,3,4,5,6], c=[10,11,12,11,13], o=['distinct', 'mean', 'min', 'collapse', 'distinct'])
        with open(results.fn, 'r') as r:
            line = r.readline()
            while line:
                sl = line.strip().split('\t')
                if len(sl[9].split(',')) == ((int(sl[2])-int(sl[1]))+1):#only keep if all bases of the motif were in the region
                    if sl[10] not in regions.keys():
                        regions[sl[10]] = [sl[0:10]]
                    else:
                        if max_motif_score_per_region:
                            if float(sl[7]) >= float(regions[sl[10]][0][7]):#get the motif that have the highest score, this can be changed to distance from the center (index 8)
                                regions[sl[10]] = [sl[0:10]]
                        else:
                            regions[sl[10]].append(sl[0:10])
                line = r.readline()
    
    with open(motifs_scored_output_file, 'w') as motifs_scored_outfile:
        for region in regions.keys():
            for motif_region in regions[region]:
                motif_region[7] = "{:.3f}".format(float(motif_region[7]))
                motifs_scored_outfile.write('\t'.join(motif_region) + '\t' + region + '\n')
        
    return regions #the highest scored motif for each tile region

def get_promoters_of_unactive_genes(genes_input_file, proms_of_unactive_genes_output_file, unactive_expr_thresh=0, gene_expression_index=-1, strand_index =5, sep='\t', num_bp_for_prom=1000):
    
    if os.path.exists(proms_of_unactive_genes_output_file):
        return proms_of_unactive_genes_output_file
    with open(genes_input_file, 'r') as genes_infile, open(proms_of_unactive_genes_output_file, 'w') as outfile:
        print 'Getting promoters (size={}) of KNOWN protein_coding genes that have expr <= {} in {}'.format(num_bp_for_prom, unactive_expr_thresh, genes_input_file)
        l = genes_infile.readline()
        while l:
            sl = l.strip().split(sep)
            if float(sl[gene_expression_index])<=unactive_expr_thresh and sl[7]=='KNOWN' and sl[8]=='protein_coding':#'pseudogene' not in sl[8]:
                if sl[strand_index] == "-":
                    outfile.write(sl[0] + sep +  sl[2] + sep + str(int(sl[2])+num_bp_for_prom) + sep + sep.join(sl[3::]) + '\n')
                else:
                    outfile.write(sl[0] + sep +  str(int(sl[1])-num_bp_for_prom)  + sep + sl[1] + sep + sep.join(sl[3::]) + '\n')
            l = genes_infile.readline()
    return proms_of_unactive_genes_output_file

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
        if so[3]=='ref' and (float(so[8])>=2 or float(so[13])>=2):
            out_file.write('\t'.join(oligos_dict[so[1]]) + '\t' + so[0] + '\t' + str(max(float(so[5]), float(so[10]))) + '\n')
    
    return output_file


def run_subset(sys_args, col_names_to_weight_param, db_name, training_dir_results, col_names):
    
    actions_list = sys_args[-1].split(',')
    prom_unactive_genes_start_index_params = 0
    if "tile_prom_regions" in actions_list:
        prom_unactive_genes_start_index_params = 7
        MPRA_tiles_input_file = sys_args[0]
        MPRA_tiles_score_per_pos_output_file= sys_args[1]
        experiment_cell_name = sys_args[2]
        cells_to_extract_info_from = sys_args[3].split(',')
        motifs_dir = sys_args[4]
        motifs_scored_output_file = sys_args[5]
        df_output_file = sys_args[6]
        if not os.path.exists(MPRA_tiles_score_per_pos_output_file):
            score_per_pos(MPRA_tiles_input_file=MPRA_tiles_input_file, output_file=MPRA_tiles_score_per_pos_output_file, experiment_cell_name=experiment_cell_name)
        if not os.path.exists(motifs_scored_output_file):
            get_motif_scores(scores_per_bp_input_file=MPRA_tiles_score_per_pos_output_file, 
                             motifs_dir=motifs_dir, motifs_scored_output_file=motifs_scored_output_file, 
                             max_motif_score_per_region=False)
        
        df_results = get_cell_info_for_motifs(motifs_scored_output_file, db_name = db_name, cells = cells_to_extract_info_from, df_output_file=df_output_file, col_names = col_names)
        df_results.to_csv(df_output_file+'.tsv', sep='\t')

    if 'prom_unactive' in actions_list:
        genes_input_file=sys_args[prom_unactive_genes_start_index_params] 
        proms_of_unactive_genes_output_file=sys_args[prom_unactive_genes_start_index_params+1]
        prom_df_output_file = sys_args[prom_unactive_genes_start_index_params+2]
        cells_to_extract_info_from_for_prom_unactive = sys_args[prom_unactive_genes_start_index_params+3].split(',')
        proms_unactive_genes = get_promoters_of_unactive_genes(genes_input_file=genes_input_file, proms_of_unactive_genes_output_file=proms_of_unactive_genes_output_file, 
                                                               unactive_expr_thresh=0, gene_expression_index=-1, strand_index =5, sep='\t', num_bp_for_prom=1000)
        prom_results_df = get_cell_info_for_regions(proms_unactive_genes, db_name = db_name, cells = cells_to_extract_info_from_for_prom_unactive, assays = ['all'], df_output_file=prom_df_output_file, col_names = col_names)
        prom_results_df.to_csv(prom_df_output_file+'.tsv', sep='\t')

    if "other_active_regions" in actions_list: 
        Other_active_regions_file = sys_args[-4]
        if not os.path.exists(Other_active_regions_file):
            Other_active_regions_file = getBed(coordinates_input_file=sys_args[-6], input_file=sys_args[-5], output_file=Other_active_regions_file)
        cells_to_extract_info_from = sys_args[-3].split(',')
        regions_df_output_file = sys_args[-2]
        Other_active_regions_df = get_cell_info_for_regions(
             Other_active_regions_file, db_name = db_name, cells = cells_to_extract_info_from, assays = ['all'], 
             sep='\t', report_cols_from_file=True, cols_indices_to_report_from_file=[4], cols_names_to_report_from_file=['Activity_Score'], df_output_file=regions_df_output_file,
             region_name_index = 3, region_strand_index = None, region_score_index = 4, motif_score_index = 4, max_number_motifs_to_report = 1,
             min_dist_from_region_start = False, min_dist_from_region_center = True, max_motif_score = False, max_region_score = False,
             col_names = col_names)

    combined_results = []    
    if "other_active_regions" in actions_list and 'tile_prom_regions' in actions_list and 'prom_unactive' in actions_list:
        combined_results =  df_results.append(prom_results_df)
        combined_results =  combined_results.append(Other_active_regions_df)    
    elif 'tile_prom_regions' in actions_list and 'prom_unactive' in actions_list:
        combined_results =  df_results.append(prom_results_df)
    elif 'other_active_regions' in actions_list and 'prom_unactive' in actions_list:
        combined_results =  Other_active_regions_df.append(prom_results_df)
    elif 'other_active_regions' in actions_list:
        combined_results =  Other_active_regions_df
    elif 'prom_unactive' in actions_list:
        combined_results =  prom_results_df
    dfout = pd.DataFrame()
    dfout_filename = "dfout_file"
    cell_name = ""
    if len(combined_results)>0:
        outcome_col = 'Activity_Score'
        cols_to_weight = []
        for c in combined_results.columns:
            if c.encode('ascii','ignore').split('___')[-1] in col_names_to_weight_param:
                cols_to_weight.append(c.encode('ascii','ignore'))
        for k in cols_to_weight:
            if '___' in k:
                dfout_filename = k.split('___')[0]+"combineddfs.tsv"
                if 'other_active_regions' in actions_list and len(actions_list)==1:
                    dfout_filename = "Other_"+k.split('___')[0]+"combineddfs.tsv"
                cell_name = k.split('___')[0]
                break
        
        dfout = combined_results.loc[:, cols_to_weight]
        dfout[outcome_col] = combined_results.loc[:, outcome_col]
        
        cell_name_series = pd.Series([cell_name for x in range(0,len(dfout))], name='CellName')
        
        dfout['CellName'] = cell_name_series#pd.concat([dfout, cell_name_series], axis=1, ignore_index=True)
        dfout_cols = {}
        for c in dfout.columns:
            if '___' in c:
                dfout_cols[c.encode('ascii','ignore')] = c.encode('ascii','ignore').split('___')[-1].lower()
            else:
                dfout_cols[c.encode('ascii','ignore')] = c.encode('ascii','ignore').lower()
        dfout.rename(columns=dfout_cols, inplace=True)
        dfout.to_csv(training_dir_results + dfout_filename, sep='\t')
    
    return dfout, dfout_filename
    
def make_binary(x, f=0):
    if x>f:
        return 1
    else:
        return 0

def make_abs(x):
    x = float(x)
    if x>=0.0:
        return x
    else:
        return x*-1.0
        
def inf_to_zero(x):

    if x == np.log2(0.0):
        return 0.0
    else:
        return x
    
def combine_lables(x):
    if 'Tss' in x or 'Tx' in x or 'BivFlnk' in x:
        return 'TSS'
    elif 'Enh' in x:
        return 'Enh'
    elif 'Repr' in x:
        return 'Repr'
    elif 'Quies' in x or 'Rpts' in x or 'Het' in x:
        return 'Quies'
    else:
        return x
        
def get_coeff(df, cols_to_weight, outcome_col, col_names_to_weight_param, dfout_filename):
    
    new_df = pd.DataFrame()
    new_df[outcome_col] = df[outcome_col].apply(make_abs).apply(make_binary, args=(0,))
    new_df[outcome_col] = new_df[outcome_col]
    for c in cols_to_weight:
        if c not in col_names_to_weight_param:
            continue
        if c=='DNase__seq'.lower() or c=='TFBinding'.lower() or c=='FANTOM'.lower():
            new_df[c] = df[c].apply(make_binary, args=(0,))
        elif c=='NumOtherTFBinding'.lower():
            new_df[c] = df[c]
        elif c=='RepliDomain'.lower() or c=='CellName'.lower() or c=='name':
            df_dummies = pd.get_dummies(df[c])
            new_df = pd.concat([new_df, df_dummies], axis=1)
        elif c == 'ChromHMM'.lower():
            df[c] = df[c].apply(combine_lables)
            df_dummies = pd.get_dummies(df[c])
            new_df = pd.concat([new_df, df_dummies], axis=1)
        elif c == 'score'.lower():
            new_df[c] = df[c]
        elif c=='TFExpr'.lower(): 
            new_df[c] = df[c].apply(float).apply(np.log2).apply(inf_to_zero)
            
    new_cols_to_weight = [c.encode('ascii','ignore') for c in new_df.columns]
    
    while 'NO' in new_cols_to_weight:#remove NO == no overlap labels 
        del new_cols_to_weight[new_cols_to_weight.index('NO')]
    
    df.to_csv(dfout_filename+'_raw.tsv', sep='\t')
    new_df.to_csv(dfout_filename, sep='\t')
    
    del new_cols_to_weight[new_cols_to_weight.index(outcome_col)]
    logit = sm.Logit(new_df[outcome_col], new_df[new_cols_to_weight]).fit(method='bfgs', maxiter=10000, full_output=True)# y=df[outcome_col], x=df[['intercept', 'HepG2___DNase__seq']])#sm.OLS(y,X).fit()#method='bfgs', full_output=True)#, maxiter=1000000000
    
    return dfout_filename, logit
    
def get_param_weights(col_names_to_weight_param, db_name, motif_info_col_names, datafiles_motifs_dir, 
                      training_dir_results, training_dir_Ernst, training_dir_Tewhey, training_dir_Vockley,
                      datafiles_HepG2_geneexpr_dir, datafiles_K562_geneexpr_dir, datafiles_GM12878_geneexpr_dir, datafiles_MCF7_geneexpr_dir):
    
    training_sets_args = ['{training_dir_Ernst}/HEPG2_SHARPR-MPRA_scores/basepredictions_HEPG2_ScaleUpDesign1and2_combinedP.txt {training_dir_results}/basepredictions_HepG2_ScaleUpDesign1and2_combinedP_perBase.txt HepG2 HepG2,Liver {datafiles_motifs_dir} {training_dir_results}/basepredictions_HepG2_ScaleUpDesign1and2_combinedP_motifs.bed {training_dir_results}/basepredictions_HepG2_ScaleUpDesign1and2_combinedP_motifs.df {datafiles_HepG2_geneexpr_dir} {training_dir_results}/HepG2_unactive_proms.bed {training_dir_results}/HepG2_unactive_proms_motifs.df HepG2,Liver tile_prom_regions,prom_unactive'.format(**locals()).split(' '),
                     '{training_dir_Ernst}/K562_SHARPR-MPRA_scores/basepredictions_K562_ScaleUpDesign1and2_combinedP.txt {training_dir_results}/basepredictions_K562_ScaleUpDesign1and2_combinedP_perBase.txt K562 K562,Whole_Blood {datafiles_motifs_dir} {training_dir_results}/basepredictions_K562_ScaleUpDesign1and2_combinedP_motifs.bed {training_dir_results}/basepredictions_K562_ScaleUpDesign1and2_combinedP_motifs.df {datafiles_K562_geneexpr_dir} {training_dir_results}/K562_unactive_proms.bed {training_dir_results}/K562_unactive_proms_motifs.df K562,Whole_Blood tile_prom_regions,prom_unactive'.format(**locals()).split(' '),
                     '{training_dir_Tewhey}/master.geuvadis.probes.oligo.sort.bed {training_dir_Tewhey}/20150102_Geuv.HepG2.txt {training_dir_results}/Tewhey_Cell2016_20150102_Geuv.HepG2_SigRegs.bed HepG2,Liver {training_dir_results}/HepG2Motifs_Tewhy other_active_regions'.format(**locals()).split(' '),
                     '{datafiles_GM12878_geneexpr_dir} {training_dir_results}/GM12878_unactive_proms.bed {training_dir_results}/GM12878_unactive_proms_motifs.df GM12878,Whole_Blood {training_dir_Tewhey}/master.geuvadis.probes.oligo.sort.bed {training_dir_Tewhey}/20150102_Geuv.12878.txt {training_dir_results}/Tewhey_Cell2016_20150102_Geuv.12878_SigRegs.bed GM12878,Whole_Blood {training_dir_results}/GM1287819239Motifs prom_unactive,other_active_regions'.format(**locals()).split(' '),
                     'None None {training_dir_Vockley}/Vockley_Cell_2016_STable2_SigRegions_75bp.bed A549,IMR-90,Lung {training_dir_results}/A549Motifs other_active_regions'.format(**locals()).split(' '),
                     '{datafiles_MCF7_geneexpr_dir} {training_dir_results}/MCF-7_unactive_proms.bed {training_dir_results}/MCF-7_unactive_proms_motifs.df MCF-7,Breast prom_unactive'.format(**locals()).split(' ')
                     ]
    files_to_use_for_training = []
    files_to_use_for_training_tsv = []
    
    for training_set_args in training_sets_args:
        df, df_tsv = run_subset(training_set_args, col_names_to_weight_param, db_name, training_dir_results+'/', motif_info_col_names[:])
        files_to_use_for_training.append(df)
        files_to_use_for_training_tsv.append(df_tsv)
    
    #print 'List of files used for training: ', ','.join(files_to_use_for_training_tsv)
    
    combined_results = files_to_use_for_training[0].append(files_to_use_for_training[1::], ignore_index=True)
    
    reported_col_names = combined_results.columns 
    outcome_col = 'Activity_Score'.lower()
    cols_to_weight = [c.encode('ascii','ignore') for c in reported_col_names]
    del cols_to_weight[cols_to_weight.index(outcome_col)]
    
    dfout_filename = "{0}/combineddfs.tsv".format(training_dir_results)
    combined_results.to_csv('{0}/combined_results.tsv'.format(training_dir_results), sep='\t')
    dfout_filename, logit_params = get_coeff(combined_results, cols_to_weight, outcome_col, col_names_to_weight_param, dfout_filename=dfout_filename)
    
    return logit_params


if __name__ == '__main__':
    
    db_name='regmotifs'
    
    datafiles_motifs_dir = '/data2/husen/ActiveMotifs/datafiles/Motifs/motifs_split_chr'
    
    datafiles_HepG2_geneexpr_dir = '/data2/husen/ActiveMotifs/datafiles/GeneExp/ENCODEGeneExpr/HepG2/HepG2.bed' 
    datafiles_K562_geneexpr_dir = '/data2/husen/ActiveMotifs/datafiles/GeneExp/ENCODEGeneExpr/K562/K562.bed'
    datafiles_GM12878_geneexpr_dir = '/data2/husen/ActiveMotifs/datafiles/GeneExp/ENCODEGeneExpr/GM12878/GM12878.bed'
    datafiles_MCF7_geneexpr_dir = '/data2/husen/ActiveMotifs/datafiles/GeneExp/ENCODEGeneExpr/MCF-7/MCF-7.bed'
    
    training_dir_results = '/data2/husen/ActiveMotifs/datafiles/TrainingSets/Weight_features_analysis'
    training_dir_Ernst = '/data2/husen/ActiveMotifs/datafiles/TrainingSets/Ernst_NatGen_2016'
    training_dir_Tewhey = '/data2/husen/ActiveMotifs/datafiles/TrainingSets/Tewhey_Cell2016' 
    training_dir_Vockley = '/data2/husen/ActiveMotifs/datafiles/TrainingSets/Vockley_Cell_2016'
    
    motif_info_col_names = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
    
    col_names_to_weight_param = ['ChromHMM'.lower(), 'DNase__seq'.lower(), 'FANTOM'.lower(), 'NumOtherTFBinding'.lower(), 'RepliDomain'.lower(), 'TFBinding'.lower(), 'TFExpr'.lower(), 'score'.lower()]#sys.argv[1].split(',')#
    
    logit_params = get_param_weights(col_names_to_weight_param, db_name, motif_info_col_names, datafiles_motifs_dir, 
                      training_dir_results, training_dir_Ernst, training_dir_Tewhey, training_dir_Vockley,
                      datafiles_HepG2_geneexpr_dir, datafiles_K562_geneexpr_dir, datafiles_GM12878_geneexpr_dir, datafiles_MCF7_geneexpr_dir)
    
    print logit_params.summary()
    print np.exp(logit_params.params)
    

