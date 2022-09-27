"""
Created on Sep 26, 2022

@author: Mark Melzer

Test functions to check section 1 in funMotifsMain.py
"""

import argparse
import unittest
import difflib
import sys

sys.path.insert(1, '../src/')
import DataProcessing
import Utilities




class TestSection1(unittest.TestCase):

    '''def testCollectingData1(self):
        """ Test with existing data directory """
        args = parse_args()
        params = Utilities.get_params(args.param_file)
        data_dir = DataProcessing.collect_all_data(params['all_chromatin_makrs_all_cells_combined_dir_path'],
                                                   params['data_tracks'])

        return'''

    def testCollectingData2(self):
        """ Test without existing data directory """
        data_tracks = "./InputTestFilesSection1/DataTracks/CAGE_expr_per_peak_all_cells_promoters_enhancers.bed4,./InputTestFilesSection1/DataTracks/RoaDomainsAllGrouped.bed4,./InputTestFilesSection1/DataTracks/RoaLoopsAllGrouped.bed4,./InputTestFilesSection1/DataTracks/ReplicationDomains.bed4,./InputTestFilesSection1/DataTracks/*ChIP-seq.bed4,./InputTestFilesSection1/DataTracks/*_DNase-seq.bed4,./InputTestFilesSection1/DataTracks/*_ChromatinStates.bed4"
        all_chromatin_makrs_all_cells_combined_dir_path2='./ InputTestFilesSection1 / chromatin_marks_all_cells_onlynarrowpeaks2'
        data_dir = DataProcessing.collect_all_data(all_chromatin_makrs_all_cells_combined_dir_path2, data_tracks)
        # check if the created file is the expected output

        with open('InputTestFilesSection1/chromatin_marks_all_cells_onlynarrowpeaks/chr10.bed', 'r') as a, open(data_dir + '/chr10.bed', 'r') as b:
            differ = difflib.Differ()
            for line in differ.compare(a.readlines(), b.readlines()):
                print(line)
                self.assertNotEqual(line[0], '-')
                self.assertNotEqual(line[0], '+')
        return


if __name__ == '__main__':
    unittest.main()