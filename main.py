#!/usr/bin/env python # -*- coding: utf-8 -*-

""" Program to identify ORF in a Fasta file for the BAM course at the TU Berlin """

import re
import numpy as np
import pandas as pd


class FileReader:
    """
    Class to get a file name, open it and prepare it for further analysis
    """
    def __init__(self, path):
        self.path = path
        self.sequence = self.delete_non_seq_text(self.read_file())

    def read_file(self):
        file = open(self.path)
        sequence = file.read()
        file.close()
        return sequence

    def delete_non_seq_text(self, raw_file: str) -> str:
        text = re.findall(r'([ATGC]{10}[ATGC]*)', raw_file)
        sequence = ''.join(text)
        return sequence


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    a = FileReader(path='./assets/sequence.fasta')
    # print(a.read_file()[0:100])
    print(a.sequence)
