#!/usr/bin/python2.7

import sys
import numpy as np

seqA = "AGTAAACCGTA"
seqB = "AGTAACGTA"

class LocalAlignment:

    def __init__(self, seq1, seq2, seq_match=3, seq_miss=-3, seq_del=-2):

        self.seq1 = seq1
        self.seq1L = list(seq1)
        self.seq2L = list(seq2)
        self.seq2 = seq2
        self.seq_match = seq_match
        self.seq_miss = seq_miss
        self.seq_del = seq_del

    def start(self):
        """Initialise local sequence alignment algorithms"""

        self.boolMat = self.boolMatrixGen()


    def baseMatrix(self):
        """Generates initial base numpy matrix"""
        return np.zeros([len(self.seq1)+1, len(self.seq2)+1])


    def boolMatrixGen(self):
        """Generates boolean similarity matrix from nucleotide sequences"""

        matrix = np.zeros([len(self.seq1), len(self.seq2)], dtype=bool)

        for r in range(len(self.seq1L)):
            for c in range(len(self.seq2L)):
                if self.seq1L[r] == self.seq2L[c]:
                    matrix[r][c] = True

        return matrix


    def scoreMatrix(self):
        """Generates the score matrix from the base matrix"""

        matrix = np.zeros([len(self.seq1)+1, len(self.seq2)+1])














x = LocalAlignment(seqA, seqB)
x.start()
print(x.boolMat)
