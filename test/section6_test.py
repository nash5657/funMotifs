"""
Created on Nov 30, 2022

@author: Mark Melzer

Test functions to check section 6 in funMotifsMain.py
"""
import unittest
import sys
import os
import difflib
from glob import glob

sys.path.insert(1, '../src/')

import FindMotifMutations


def compare_files(fileA, fileB):
    """
    Function to compare an output file with the expected output to assert its correctness.
    """
    with open(fileA, 'r') as a, open(fileB, 'r') as b:
        differ = difflib.Differ()
        for line in differ.compare(a.readlines(), b.readlines()):
            if line[0] in ['-', '+']:
                print(line)
                return False
    return True


class TestSection6(unittest.TestCase):
    """
    Class to test the different functions in the file FindMotifMutations.py
    The necessary data base table can be created using section6_postgres.csv
    """

    def test_get_BedTool_from_variant_file(self):
        """
        Function to test the functionality of get_BedTool_from_variant_file
        """
        FindMotifMutations.get_BedTool_from_variant_file("InputTestFilesSection6/Variant_test_input_file.maf")
        os.remove("InputTestFilesSection6/Variant_test_input_file.maf_tmp")

        assert compare_files("InputTestFilesSection6/Variant_test_input_file.maf.bed",
                             "InputTestFilesSection6/Variant_test_expected_file.bed")

        return

    
    def test_get_BedTool_for_functional_motif(self):
        """
        Function to test the functionality of get_BedTool_for_functional_motifs
        This simultaneously tests get_BedTool_from_dataframe
        """
        funMotifs = {'liver': [1, 3]}
        tissue = "liver"
        db_user_name = "mm99"
        db_name = "funmotifsdb"
        output_file = "InputTestFilesSection6/funMotif_test_output_file.bed"

        FindMotifMutations.get_BedTool_for_functional_motifs(funMotifs, tissue, db_user_name, db_name, output_file)

        assert compare_files("InputTestFilesSection6/funMotif_test_output_file.bed",
                             "InputTestFilesSection6/funMotif_test_expected_file.bed")
        return

    def test_find_funMotif_variants_in_tissue(self):
        """
        Function to test the functionality of find_funMotif_variants_in_tissue
        This simultaneously test overlap_variants_and_motifs
        """
        variant_file = "InputTestFilesSection6/Variant_test_input_file.maf"
        FindMotifMutations.get_BedTool_from_variant_file(variant_file)
        funMotifs = {'liver': [1, 2, 3]}
        tissue = 'liver'
        db_name = "funmotifsdb"
        db_user_name = "mm99"
        output_file = "InputTestFilesSection6/Overlap_Variants_Motifs"
        motif_BedTool_file = "InputTestFilesSection6/motif_BedTool_file.bed"
        FindMotifMutations.find_funMotif_variants_in_tissue(funMotifs, tissue, variant_file + ".bed", db_name,
                                                            db_user_name, output_file, motif_BedTool_file)

        assert compare_files("InputTestFilesSection6/Overlap_Variants_Motifs",
                             "InputTestFilesSection6/Overlap_Variants_Motifs_expected")

        return


if __name__ == '__main__':
    prev = glob("InputTestFilesSection6/*")
    unittest.main(exit=False)

    rm = []
    after = glob("InputTestFilesSection6/*")
    for i in after:
        if i not in prev:
            rm.append(i)
    for i in rm:
        os.remove(i)
