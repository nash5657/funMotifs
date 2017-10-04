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
    



