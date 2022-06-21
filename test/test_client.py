#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import re
from main import FileReader


class TestClient(unittest.TestCase):
    """
    Class that contains all the test necessary to make the program
    """

    dummies = {'text_1': '> NZ_CP027599.1Escherichia coli strain 97 - 3250 chromosome, complete genome ATCCCGGCCCCGG' \
                       'CAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTG\nAAGGTAAATCTAACCAACTGGCGCGCGCGGCGG' \
                       'CTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAA',
               'text_2': 'CCCCCATGTGGACCCCATAGCCCCCCATGGCGAGACCGTGTATGTAACCCCCATGATGATGATGATGTGACCCCTTGTGACCCCCCCGTGTTT' \
                         'GCCTAACCCCCCATAGACGACCGACGAAGATGACCCCCATGAGACGCGAGCGAGCGTAA'}

    def test_get_file_name_return_str(self):
        """
        Test to prove that the function read_file() form class FileReader can read the file name and return the content
        as a str.
        """
        ner = FileReader('../test/dummy.fasta')
        reader = ner.read_file()
        self.assertIsInstance(reader, str)

    def test_delete_non_sequence_text(self):
        nero = FileReader('../test/dummy.fasta')
        # non sequence text plus additional 2 strings
        ner = nero.delete_non_seq_text(TestClient.dummies['text_1'])
        # Only sequence text with multiple new line characters
        ner2 = nero.sequence

        self.assertEqual('ATCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAAC' \
                         'TGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAA', ner)
        self.assertEqual('ATCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTG' \
                         'AAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAA' \
                         'CCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTAACGGCATT' \
                         'ATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACATGGTTAAAG' \
                         'CCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTGATCGACGA' \
                         'TATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCATACCTTCAACGCCCTGCTGGAA' \
                         'GGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTTTGA' \
                         'AATCCCGCTTCGGTTGGGGACTGACTGTGGCAATCGAACCGCCAGAGCTGGAAACCCGTGTGGCGATCCT' \
                         'GATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCAAGCGTCTA' \
                         'CGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTACCGGACGGG' \
                         'CGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTGGTCACCAT' \
                         'CGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCAAGCGTCGA' \
                         'TCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGC', ner2)

    def test_orf_identifier(self):
        ner = FileReader('../test/dummy.fasta')
        orf_list = ner.orf_finder(TestClient.dummies['text_2'], 'eukaryote')
        self.assertEqual(['ATGTGGACCCCATAG', 'ATGGCGAGACCGTGTATGTAA', 'ATGATGATGATGATGTGA',
                          'ATGACCCCCATGAGACGCGAGCGAGCGTAA'],
                         orf_list)
        # self.assertEqual(['ATGTGGACCCCATAG', 'ATGGCGAGACCGTGTATGTAA', 'ATGATGATGATGATGTGA', 'TTGTGA', 'GTGTTTGCCTAA',
        #                   'ATAGACGACCGACGAAGATGA', 'ATGAGACGCGAGCGAGCGTAA'], orf_list)

