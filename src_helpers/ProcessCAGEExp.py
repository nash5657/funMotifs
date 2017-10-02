'''
Created on May 4, 2016

@author: Husen M. Umer
'''
import sys
import numpy as np
'''
output format: peak1 info(chr start end strand) cell1#avgExpValue,cell2#avgExpValue,celln#avgExpValue
'''
def extract_expression_per_peak_per_cell(cell_to_sampleID_mapping_inputfile, 
                                         peak_expression_matrix_inputfile, 
                                         start_index_expr_values=7, 
                                         peak_coords_index = 0, 
                                         header_line_starting_keyword="00Annotation", 
                                         peak_lines_starting_keyword="chr", 
                                         ignore_lines_starting_wiht_keyword='#'):
    #peak coord format in FANTOM dataset =  chr10:100013403..100013414,-
    cell_sample_dict = {}
    with open(cell_to_sampleID_mapping_inputfile) as cell_to_sampleID_mapping_infile:
        line = cell_to_sampleID_mapping_infile.readline()
        while line!="":
            cell_sample_dict[line.strip().split('\t')[0]] = line.strip().split('\t')[1].split(',')
            line = cell_to_sampleID_mapping_infile.readline()
    
    peak_expression_matrix_outputfile = peak_expression_matrix_inputfile+"_avgExprValueperCell.bed4"
    peak_expression_matrix_outfile = open(peak_expression_matrix_outputfile, 'w')
    
    with open(peak_expression_matrix_inputfile, 'r') as peak_expression_matrix_infile:
        line = peak_expression_matrix_infile.readline()
        sample_ids_from_expr_file = []
        while line != "":
            if not line.strip().startswith(ignore_lines_starting_wiht_keyword):
                split_line = line.split('\t')
                if line.startswith(header_line_starting_keyword) and len(sample_ids_from_expr_file) == 0: #the header line starts with keyword
                    for x in split_line[start_index_expr_values::]:
                        sample_ids_from_expr_file.append(x.strip().split('.')[2]) #an example of a sample ID is: tpm.293SLAM%20rinderpest%20infection%2c%2000hr%2c%20biol_rep3.CNhs14408.13543-145H6 but since only the libraryID (e.g CNhs14408) is listed in the sample_info file so only that one is extracted 

                elif line.startswith(peak_lines_starting_keyword):
                    new_line = []
                    for cell in cell_sample_dict.keys():
                        list_of_expression_in_this_peak_from_all_samples_of_this_cell = []
                        for sampleID in cell_sample_dict[cell]:
                            list_of_expression_in_this_peak_from_all_samples_of_this_cell.append(float(split_line[start_index_expr_values+sample_ids_from_expr_file.index(sampleID)]))
                        #only add cell lines that have the peak active in at least one of their samples; Note: changes the next line to get only cells that are active in a given percentage of samples because the list contains the value from all samples of this cell
                        if np.mean(list_of_expression_in_this_peak_from_all_samples_of_this_cell) > 0: #=0 means write all regardless of their activity
                            new_line.append(cell + "#FANTOM#" + str(np.mean(list_of_expression_in_this_peak_from_all_samples_of_this_cell)))
                    #only write peaks that are active in at least one cell line     
                    if len(new_line)>0:#=0 means write all the peaks regardless their activity
                        peak_expression_matrix_outfile.write(split_line[peak_coords_index].strip().split(',')[0].replace(':',"\t").replace('..','\t') + '\t' + ','.join(new_line) + '\n')
            line = peak_expression_matrix_infile.readline()
        peak_expression_matrix_outfile.close()
    return peak_expression_matrix_outputfile

if __name__ == '__main__':
    cell_to_sampleID_mapping_inputfile = sys.argv[1]
    peak_expression_matrix_inputfile = sys.argv[2]
    extract_expression_per_peak_per_cell(cell_to_sampleID_mapping_inputfile, peak_expression_matrix_inputfile)
    