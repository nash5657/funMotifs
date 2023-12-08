import re

def extract_variables_from_file(file_path):
    """ Extracts variables and their values from a file and returns them as a dictionary. """
    variable_dict = {}

    # Regular expression pattern to identify variable names and values
    pattern = r'([\w/,]+)=([-+]?\d*\.\d+|\d+)'

    with open(file_path, 'r') as file:
        for line in file:
            # Find all matches in the line
            matches = re.findall(pattern, line)
            for vars, val in matches:
                # Splitting multiple variables if they exist
                vars_split = vars.split(',')
                for var in vars_split:
                    # Adding to the dictionary
                    variable_dict[var] = float(val)

    return variable_dict

# Usage
file_path = 'annotation_wights.txt'  # Replace with the path to your file
variables = extract_variables_from_file(file_path)
values_list = list(variables.values())
