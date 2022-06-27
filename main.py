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

    def orf_finder(self, sequence: str, organism: str) -> list:
        """
        Function that use regex to find the posible ORFs in a sequence. The function can take different sequences and
        organism types as prokaryote, eukaryote and mitochondria that have specific start codons. The function scans the
        sequence one base at a time until it finds a start codon and the findall method of the regex library takes non-
        overlapping sequence matches. The matches are stored in a list and retured
        :sequence type=str Sequence to be analyzed
        :organism type=str Specific organism that corresponds to the sequence. (It affects the start coddons)
        :return: type=list  List of ORFs
        """
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

        # https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/start-codon
        c = 0
        tmp_orf = []
        loop = True
        append_var = False

        while loop:
            if sequence[c + 0] + sequence[c + 1] + sequence[c + 2] in FileReader.start_codon[organism] and append_var == False:
                append_var = True
                tmp_orf.append('ATG')
                c += 3
            elif append_var and sequence[c + 0] + sequence[c + 1] + sequence[c + 2] in {'TAG', 'TAA', 'TGA'}:
                append_var = False
                tmp_orf.append(sequence[c + 0] + sequence[c + 1] + sequence[c + 2]+',')
                c += 3
            elif append_var:
                tmp_orf.append(sequence[c + 0] + sequence[c + 1] + sequence[c + 2])
                c += 3
            else:
                c += 1


            if c+2 >= len(sequence):
                loop = False

        orf_list = ''.join(tmp_orf)
        orf_list = orf_list.split(',')
        try:
            orf_list.remove('')
        except:   # If it doesn't have the empty '' it means that the segment doesn't have a termination codon
            last_codons = FileReader.orf_finder_slow(self, orf_list[-1][1:], organism)
            orf_list.remove(orf_list[-1])
            for i in last_codons:
                orf_list.append(i)
            print(last_codons)

        return orf_list


if __name__ == '__main__':
    a = FileReader(path='./assets/sequence.fasta')
    # print(a.read_file()[0:100])
    print(a.orf_finder_slow(a.sequence, 'eukaryote'))
