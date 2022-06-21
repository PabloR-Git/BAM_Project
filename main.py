#!/usr/bin/env python # -*- coding: utf-8 -*-

""" Program to identify ORF in a Fasta file for the BAM course at the TU Berlin """

import re
import numpy as np
import pandas as pd


class FileReader:
    """
    Class to get a file name, open it and prepare it for further analysis
    """
    start_codon = {'eukaryote': {'ATG'},
                   'prokaryote': {'ATG', 'GTG', 'TTG'},
                   'mitochondria': {'ATG', 'ATA'}}

    def __init__(self, path):
        self.path = path
        self.sequence = self.delete_non_seq_text(self.read_file())

    def read_file(self):
        """
        RAW File reader function.
        """
        file = open(self.path)
        sequence = file.read()
        file.close()
        return sequence

    def delete_non_seq_text(self, raw_file: str) -> str:
        """
        Function to delete the non sequence characters contained in the RAW files.
        :raw_file type=str Raw text without non sequence header
        :return: string without new line spaces ready to analyze
        """
        text = re.findall(r'([ATGC]{10}[ATGC]*)', raw_file)
        sequence = ''.join(text)
        return sequence

    def orf_finder(self, sequence, organism):
        if organism == 'eukaryote':
            orf_list = re.findall(r'((ATG)(([ATGC]{3})*?)(TAG|TAA|TGA))', sequence)
        if organism == 'prokaryote':
            orf_list = re.findall(r'((ATG|GTG|TTG)(([ATGC]{3})*?)(TAG|TAA|TGA))', sequence)
        if organism == 'mitochondria':
            orf_list = re.findall(r'((ATG|ATA)(([ATGC]{3})*?)(TAG|TAA|TGA))', sequence)

        for i in range(len(orf_list)):
            orf_list[i] = orf_list[i][0]
        return orf_list

    def orf_finder_slow(self, sequence: str, organism: str) -> list:
        # TODO: Change the start codon find strategy from 3 to one !!!!
        # https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/start-codon
        c = 0
        tmp_orf = []
        loop = True
        append_var = False

        while loop:
            if sequence[c + 0] + sequence[c + 1] + sequence[c + 2] in FileReader.start_codon[organism] and append_var == False:
                append_var = True
                tmp_orf.append('ATG')
            elif sequence[c + 0] + sequence[c + 1] + sequence[c + 2] in {'TAG', 'TAA', 'TGA'} and append_var == True:
                append_var = False
                tmp_orf.append(sequence[c + 0] + sequence[c + 1] + sequence[c + 2]+',')
            elif append_var == True:
                tmp_orf.append(sequence[c + 0] + sequence[c + 1] + sequence[c + 2])
            else:
                pass

            c += 3

            if c+2 >= len(sequence):
                loop = False

        orf_list = ''.join(tmp_orf)
        orf_list = orf_list.split(',')
        orf_list.remove('')

        return orf_list


if __name__ == '__main__':
    a = FileReader(path='./assets/sequence.fasta')
    # print(a.read_file()[0:100])
    print(a.orf_finder_slow(a.sequence, 'eukaryote'))
