'''
Created on 28 Sep 2017

@author: husensofteng
'''
import json, os

def retreive_TFFamilyName_for_motifNames(TF_family_matches_file):#TFFamilyName TF_name
    "Retrieves the TF family name for each TF name"
    
    motifTFName_TFNames_matches_dict = {}
    with open(TF_family_matches_file, 'r') as TFFamily_matches_infile:
        lines = TFFamily_matches_infile.readlines()
        for line in lines:
            sl = line.strip().split('\t')
            if len(sl)>1:
                motifTF_name = sl[0].strip().upper()
                if motifTF_name not in motifTFName_TFNames_matches_dict.keys():
                    motifTFName_TFNames_matches_dict[motifTF_name] = []
                for s in sl:
                    if s.strip().upper()!="": 
                        if s.strip().upper() not in motifTFName_TFNames_matches_dict[motifTF_name]:
                            motifTFName_TFNames_matches_dict[motifTF_name].append(s.strip().upper())                       
    return motifTFName_TFNames_matches_dict

#Given a GTEX file retrieve gene expression from each tissue for each TF name (as defined in the retreive_TFFamilyName_for_motifNames) 

def get_expression_level_per_originType_per_TF(motifTFName_TFNames_matches_dict, 
                                               normal_gene_expression_inputfile,
                                               origin_gene_expression_values_outputfile, 
                                               index_tissues_names_row_start, 
                                               index_gene_names_col,
                                               index_gene_values_start, 
                                               sep):
    tissue_origin_gene_expression_values = {}
    if not os.path.exists(origin_gene_expression_values_outputfile):
        tf_names_to_extract_gene_expression_for = []#list_tf_names_from_tracks#get names of TFs from the TFFamily file and the dirs contaning ChIP-seq datasets
        for k in motifTFName_TFNames_matches_dict.keys():
            tf_names_to_extract_gene_expression_for.append(k)
        tf_names_to_extract_gene_expression_for = list(set(tf_names_to_extract_gene_expression_for))
        TFs_extract_expression = get_TFs_extract_expression(tf_names_to_extract_gene_expression_for)
        with open(normal_gene_expression_inputfile) as normal_gene_expression_infile:
            for i in range(0, index_tissues_names_row_start):
                line = normal_gene_expression_infile.readline()#skip until it gets to the tissue names row
            tissue_names = normal_gene_expression_infile.readline().strip().split(sep)[index_gene_names_col+1::]
            line = normal_gene_expression_infile.readline()
            
            #iniatilze the tissue_origin_gene_expression_values matrix
            for tissue in tissue_names:
                tissue_origin_gene_expression_values[tissue] = {}
                for tf in TFs_extract_expression:
                    tissue_origin_gene_expression_values[tissue][tf] = 'NaN'
            while line:
                gene_name = line.strip().split(sep)[index_gene_names_col].upper()
                gene_values = line.strip().split(sep)[index_gene_names_col+1::]
                motif_tf_names_corresponding_to_this_gene = get_tfName_fromGeneName(TFs_extract_expression, gene_name)
                if len(motif_tf_names_corresponding_to_this_gene)>0:
                    for tf_name in motif_tf_names_corresponding_to_this_gene:
                        for i in range(0, len(tissue_names)):
                            if tissue_origin_gene_expression_values[tissue][tf_name] == 'NaN':
                                tissue_origin_gene_expression_values[tissue_names[i]][tf_name] = float(gene_values[i])
                            else:
                                tissue_origin_gene_expression_values[tissue_names[i]][tf_name] += float(gene_values[i])
                line = normal_gene_expression_infile.readline()
        #write the results to output file
        with open(origin_gene_expression_values_outputfile, 'w') as origin_gene_expression_values_outfile:
            json.dump(tissue_origin_gene_expression_values,origin_gene_expression_values_outfile)
            #for tissue_name in tissue_origin_gene_expression_values:
            #    for gene_name in tissue_origin_gene_expression_values[tissue_name]:
            #        origin_gene_expression_values_outfile.write(tissue_name + sep + gene_name + sep + str(tissue_origin_gene_expression_values[tissue_name][gene_name]) + '\n')
    else:
        with open(origin_gene_expression_values_outputfile, 'r') as origin_gene_expression_values_infile:
            tissue_origin_gene_expression_values = json.load(origin_gene_expression_values_infile)
    return tissue_origin_gene_expression_values

def get_tfName_fromGeneName(TFs_extract_expression, gene_name):
    motif_tf_names_corresponding_to_this_gene = []
    for tf in TFs_extract_expression:
        if gene_name==tf or gene_name in TFs_extract_expression[tf]:
            motif_tf_names_corresponding_to_this_gene.append(tf)
    return motif_tf_names_corresponding_to_this_gene

def get_TFs_extract_expression(list_of_tf_names):
    TFs_extract_expression = {}
    for tf_name in list_of_tf_names:
        if tf_name not in TFs_extract_expression:
            TFs_extract_expression[tf_name] = []
            if '(' in tf_name:
                TFs_extract_expression[tf_name].append(tf_name.split('(')[0])
            #add TF names after removing dashes and splitting them by :: (in case they were in the name)
            if '-' in tf_name:#in case - was in the name remove it
                TFs_extract_expression[tf_name].append(tf_name.replace('-',''))
            if '::' in tf_name:
                s_tf_name = tf_name.split('::')
                for s in s_tf_name:
                    if s not in TFs_extract_expression[tf_name]:
                        TFs_extract_expression[tf_name].append(s)
    return TFs_extract_expression
