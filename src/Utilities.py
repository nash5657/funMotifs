'''
Created on 28 Sep 2017

@author: husensofteng
@contributor: Mark Melzer
Utility functions
-set a temp dir for bedtools intermediate files
-get parameters from the sys.argv and the given conf file
'''
import pandas as pd


def get_value(str):
    if str in ["yes", "YES", "Yes", "y", "Y", "True", "true", "TRUE", "1", 1, True]:
        return True
    else:
        return False


def get_params(param_file):
    params = {}

    print(param_file)
    with open(param_file, 'r') as params_infile:
        params_from_file = params_infile.readlines()
        for line in params_from_file:
            if not line.startswith('//') and not line.startswith('#') and '=' in line:
                if len(line.strip().split('=')) == 2:
                    if line.strip().split('=')[0] not in list(params.keys()):
                        params[line.strip().split('=')[0]] = line.strip().split('=')[1].split('#')[0].split('//')[0]
    return params


def retreive_key_values_from_dict_file(dict_input_file, key_value_sep, values_sep):  # TFFamilyName TF_name
    """
    Retrieves the key and its values"
    """
    # list to save representative cell names
    rep_cell_names = []
    # dictionary to map cell name onto representative cell name
    value_key_dict = {}
    with open(dict_input_file, 'r') as dict_input_file_infile:
        lines = dict_input_file_infile.readlines()
        for line in lines:
            if line.startswith('//') or line.startswith('#'):  # or '=' not in line:
                continue
            if not line.__contains__(key_value_sep):
                rep_cell_names.append(line)
            else:
                sl = line.strip().split(key_value_sep)
                rep_cell_names.append(sl[0].strip())
                if not sl[1].__contains__(values_sep):
                    if sl[1].strip() not in list(value_key_dict.keys()):
                        value_key_dict[sl[1].strip()] = sl[0].strip()
                else:
                    values = sl[1].strip().split(values_sep)
                    for value in values:
                        if value not in list(value_key_dict.keys()):
                            value_key_dict[value] = sl[0].strip()

    return rep_cell_names, value_key_dict


def cell_to_tissue_matches(dict_input_file, key_value_sep='=', value_sep=','):
    # TODO: write unit test
    value_key_dict = {}
    with open(dict_input_file, 'r') as infile:
        lines = infile.readlines()
        for line in lines:
            if not line.__contains__(key_value_sep):
                continue
            sl = line.split(key_value_sep)
            values = sl[1].split(value_sep)
            for value in values:
                value_key_dict[value.strip('\n')] = sl[0]
    return value_key_dict


def bed_to_cell_wise_dataframe(bed_file, cells_to_assay, output_columns=None):
    # TODO: write unit test
    """
    Function takes bed file and transforms it to data frame.
    Motifs that appear in several cells will be in a separate line per cell
    """
    if output_columns is None:
        output_columns = ['posrange', 'cell', 'chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
    df = pd.read_csv(bed_file, sep='\t')
    cells = []
    for c in df.columns:
        if c not in output_columns:
            if not c.__contains__('___'):
                output_columns.append(c)
            else:
                tmp = c.split('___')
                cells.append(tmp[0])
                assay = tmp[1]
                if assay not in output_columns:
                    output_columns.append(assay)
    new_df = pd.DataFrame(columns=output_columns)
    for ind in df.index:
        for cell in cells:
            tmp_df = {'cell': cell}
            for col in output_columns:
                if col in df.columns:
                    tmp_df[col] = df[col][ind]
            # TODO: check for upper or lowercase letters and prevent errors
            for assay in cells_to_assay[cell]:
                tmp_df[assay] = df[cell + '___' + assay][ind]
        new_df.append(tmp_df)
    new_df.to_csv('...')
    return new_df


'''
desc: reports motif-breaking info for mutations in each motif sites. It check the difference between nucleotide frequency of the ref allele and the mutated-to allele of each mutation in the PWM of the anchor motif.
in: motifs_PFM matrix (get it from ENCODE, for instance), mutated motifs infput file (contains mutations info), output file name
out: all mutated motifs with adding breaking info to the mutations, all mutations at the motifs with adding breaking info  
out: only mutated motifs have motif-breaking mutation (difference between ref and mut allele > given_Threshold) and the fourth returned file is the list of motif breaking mutations (diff>Threshold)
'''


def get_freq_per_motif(motif_PFM_input_file):
    "given a PFM file return a dict, a key for each tf and the freq as value"
    PFM_motifs_lines = [""]
    with open(motif_PFM_input_file, 'r') as PFM_motifs_infile:
        PFM_motifs_lines = PFM_motifs_infile.readlines()

    PFM_motifs_dict = {}
    nucleotides = ['A', 'C', 'G', 'T']  # default nucleotides
    motif_info_sep = ' '
    motif_name = ""
    if 'MEME' in PFM_motifs_lines[0]:
        motif_info_sep = ' '
        for line in PFM_motifs_lines:
            if 'ALPHABET=' in line:
                nucleotides = []
                ALPHABETS = line.split('=')[1].strip()
                for alph in ALPHABETS:
                    nucleotides.append(alph.upper())

            if line.strip() != "" and not line.startswith('letter') and not line.startswith('URL'):
                if line.startswith('MOTIF'):
                    motif_name = line.strip().split(motif_info_sep)[2] + '_' + line.strip().split(motif_info_sep)[1]
                else:
                    if motif_name != "":  # if it has been initialized
                        if motif_name not in PFM_motifs_dict:
                            PFM_motifs_dict[motif_name] = []
                        split_line = line.strip().split()
                        freq_per_allele = []
                        for s in split_line:
                            try:
                                freq_per_allele.append(float(s.strip()))
                            except ValueError:
                                continue
                        if len(freq_per_allele) == len(nucleotides):  # freq of the 4 nucleotides
                            nucl_weigts = {}
                            for i, nucl in enumerate(nucleotides):
                                nucl_weigts[nucl] = float(freq_per_allele[i])
                            PFM_motifs_dict[motif_name].append(
                                nucl_weigts)  # , C: float(split_line[1]), G: float(split_line[2]), T: float(split_line[3])})
    return PFM_motifs_dict
