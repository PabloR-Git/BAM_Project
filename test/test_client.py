#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from main import FileReader


class TestClient(unittest.TestCase):
    """
    Class that contains all the test necessary to make the program
    """

    dummies = {'text_1': '> NZ_CP027599.1Escherichia coli strain 97 - 3250 chromosome, complete genome ATCCCGGCCCCGG' \
                       'CAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTG\nAAGGTAAATCTAACCAACTGGCGCGCGCGGCGG' \
                       'CTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAA'}

    def test_get_file_name_return_str(self):
        ner = FileReader('../assets/dummy.fasta')
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
