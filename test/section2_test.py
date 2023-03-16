"""
Created on Sep 28, 2022

@author: Mark Melzer

Test functions to check section 2 in funMotifsMain.py
"""

import unittest
import sys

sys.path.insert(1, '../src/')

import MotifAnnotation


class TestSection2(unittest.TestCase):

    def testFunction(self):
        normal_expression = {u'Thyroid': {u'EGR1': 22.64, u'ZFX': 0.1732, u'BACH1::MAFK': 4.558}}
        matching = {'Thyroid': ['Thyroid'], 'HepG2': ['HepG2'], 'GSS': ['GSS'], 'GM12878': ['GM12878']}
        motifTFName = {'NR1H3::RXRA': ['NR1H3::RXRA', '#MG', 'RXRA'], 'OTX1': ['OTX1', '#M'], 'FOXQ1': ['FOXQ1', '#M'],
                       'GSC2': ['GSC2', '#M'], 'SP1': ['SP1', '#ME'], 'SP2': ['SP2', '#ME']}
        cell_assays = {u'Thyroid': {u'TFExpr': 0.0}, u'HepG2': {u'RepliDomain': 'NO', u'FANTOM': 0.0},
                       u'GSS': {u'FANTOM': 0.0}, u'GM12878': {u'ContactingDomain': 0.0, u'LoopDomain': 0.0}}
        cell_tfs = {u'GM12878': [u'ZNF143', u'SP1', u'RFX5', u'MAZ', u'MTA3', u'BATF', u'POLR2AphosphoS5', u'EED']}
        tf_cells = {u'MAZ': [u'GM12878'], u'POLR2AphosphoS5': [u'GM12878'], u'EED': [u'GM12878'],
                    u'ZNF143': [u'GM12878'],
                    u'MTA3': [u'GM12878'], u'SP1': [u'GM12878'], u'RFX5': [u'GM12878'], u'BATF': [u'GM12878']}
        assay_cells = {u'TFExpr': u'real', u'NumOtherTFBinding': u'real', u'ContactingDomain': u'real',
                       u'OtherTFBinding': u'text', u'LoopDomain': u'real', u'DNase-seq': u'real', u'ChromHMM': u'text',
                       u'TFBinding': u'real', u'RepliDomain': u'text', u'FANTOM': u'real'}
        motifs_overlapping_tracks_files, scored_motifs_overlapping_tracks_files = MotifAnnotation.run_overlay_resources_score_motifs(
            './InputTestFilesSection2/motifs_per_chr',
            './InputTestFilesSection2/chromatin_marks_all_cells_onlynarrowpeaks',
            './InputTestFilesSection2/overlapping_scored_motifs_onlynarrowpeaks',
            'yes',
            12,
            normal_expression,  # output file
            matching,  # output
            motifTFName,  # output
            cell_assays,  # output
            cell_tfs,  # output
            tf_cells,  # output
            assay_cells)  # output
        print(motifs_overlapping_tracks_files, scored_motifs_overlapping_tracks_files)
        # TODO: write assertion file, so far only compared in- and output manually
        return



if __name__ == '__main__':
    unittest.main()
