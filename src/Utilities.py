'''
Created on 28 Sep 2017

@author: husensofteng
Utility functions
-set a temp dir for bedtools intermediate files
-get parameters from the sys.argv and the given conf file
'''

def get_value(str):
    if 'true' in str.lower() or 'yes' in str.lower():
        return True
    else:
        return False
    

def get_params(params_list):
    params = {}
    for arg in params_list:#priority is for the command line
        if '=' in arg: 
            if len(arg.strip().split('='))==2:
                if arg.split('=')[0] not in params.keys():
                    params[arg.strip().split('=')[0]] = arg.strip().split('=')[1]
    if 'param_file' in params:     
        with open(params['param_file'], 'r') as params_infile:
            params_from_file = params_infile.readlines()
            for line in params_from_file:
                if not line.startswith('//') and not line.startswith('#') and '=' in line:
                    if len(line.strip().split('='))==2:
                        if line.strip().split('=')[0] not in params.keys():
                            params[line.strip().split('=')[0]] = line.strip().split('=')[1].split('#')[0].split('//')[0]
    return params


def retreive_key_values_from_dict_file(dict_input_file, key_value_sep, values_sep):#TFFamilyName TF_name
    "Retrieves the key and its values"
    key_values_dict = {}
    value_key_dict = {}
    with open(dict_input_file, 'r') as dict_input_file_infile:
        lines = dict_input_file_infile.readlines()
        for line in lines:
            if line.startswith('//') or line.startswith('#'):# or '=' not in line:
                continue
            sl = line.strip().split(key_value_sep)
            key_value = sl[0].strip()
            if key_value not in key_values_dict.keys():
                key_values_dict[key_value] = []
            if key_value not in value_key_dict:
                value_key_dict[key_value]=[key_value]
            if len(sl)>1:
                for s in sl[1].split(values_sep):
                    if s.strip()!="" and s.strip() not in key_values_dict[key_value]:
                        key_values_dict[key_value].append(s.strip())
                    if s.strip()!="":
                        if s.strip() not in value_key_dict:
                            value_key_dict[s.strip()]=[]
                        value_key_dict[s.strip()].append(key_value)
    return key_values_dict, value_key_dict
    

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
    nucleotides  = ['A', 'C', 'G', 'T']#default nucleotides
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
            
            if line.strip()!="" and not line.startswith('letter') and not line.startswith('URL'):
                if line.startswith('MOTIF'):
                    motif_name = line.strip().split(motif_info_sep)[2]+'_'+line.strip().split(motif_info_sep)[1]
                else:
                    if motif_name!="":#if it has been initialized
                        if motif_name not in PFM_motifs_dict:
                            PFM_motifs_dict[motif_name] = []
                        split_line = line.strip().split()
                        freq_per_allele = []
                        for s in split_line:
                            try:
                                freq_per_allele.append(float(s.strip()))
                            except ValueError:
                                continue
                        if len(freq_per_allele)==len(nucleotides): #freq of the 4 nucleotides
                            nucl_weigts = {}
                            for i,nucl in enumerate(nucleotides):
                                nucl_weigts[nucl] = float(freq_per_allele[i])
                            PFM_motifs_dict[motif_name].append(nucl_weigts)#, C: float(split_line[1]), G: float(split_line[2]), T: float(split_line[3])})
    return PFM_motifs_dict
