#!/usr/bin/python2.7

# Author: Jonathan Hull (jonathan.hull11@gmail.com)

import sys
import random
import numpy as np

seqA = "AGTAAACCGTA"
seqB = "AGTAACGTA"

class LocalAlignment:

    def __init__(self,
                seq1, 
                seq2,
                score_method="affine",
                seq_opening=5,
                seq_gap=1,
                #seq_gap=2,
                seq_match=3, 
                seq_miss=3):

        self.seq1 = seq1
        self.seq2 = seq2
        self.score_method = score_method.lower()
        self.seq_match = seq_match
        self.seq_miss = seq_miss
        self.seq_gap = seq_gap
        self.seq_opening = seq_opening

    def build(self):
        pass


    def start(self):
        """Initialise local sequence alignment algorithms"""

        self.boolMat = self.boolMatrixGen()

        if self.score_method == "affine":
            self.scoreMat = self.AffineScoreMatrix(self.boolMat)
        elif self.score_method == "linear":
            self.scoreMat = self.scoreMatrix(self.boolMat)
        else:
            raise ComponentNameError('score_method argument must be either '
                    '"affine, or "linear".')


        print(self.scoreMat)
        pass


        pathTup = self._pathMatrix(self.scoreMat)
        base = self._seqAssembly(pathTup)
        self._seqStatistics(base)



    def boolMatrixGen(self):
        """Generates boolean similarity matrix from nucleotide sequences"""

        matrix = np.zeros([len(self.seq1), len(self.seq2)], dtype=bool)

        seq1L = list(self.seq1)
        seq2L = list(self.seq2)

        for r in range(len(seq1L)):
            for c in range(len(seq2L)):
                if seq1L[r] == seq2L[c]:
                    matrix[r][c] = True

        return matrix


    def scoreMatrix(self, boolMat):
        """Generates the score matrix from the base matrix"""
        ## Rename linear score matrix.

        matrix = np.zeros([len(boolMat)+1, len(boolMat[0])+1])

        for r in range(1, len(matrix)):
            for c in range(1, len(matrix[0])):

                if boolMat[r-1][c-1] == True:
                    diag = matrix[r-1][c-1] + self.seq_match
                else:
                    diag = matrix[r-1][c-1] - self.seq_miss

                vertical = matrix[r-1][c] - self.seq_gap
                horizontal = matrix[r][c-1] - self.seq_gap

                matrix[r][c] = sorted([diag, vertical, horizontal, 0],
                        reverse=True)[0]

        return matrix



    def AffineScoreMatrix(self, boolMat):
        """Generates scoring matrix utilising affine transformations."""

        ## Uses the equation W = u(k-1) + V.
        ## where u = gap extention penalty, V = gap opening penalty.

        ## Note: could use one list for kvc/kv and delete entries if they don't
        ## stack.

        matrix = np.zeros([len(boolMat)+1, len(boolMat[0])+1])

        kv = [] # List of gaps passed to next iteration.
        kvc = [] # tuple list :: (column index, gap extention(k))
        kvv = 1 # Vertical gap extention penalty value.
        kh = 1

        for r in range(1, len(matrix)):
            for c in range(1, len(matrix[0])):

                if boolMat[r-1][c-1] == True:
                    diag = matrix[r-1][c-1] + self.seq_match
                else:
                    diag = matrix[r-1][c-1] - self.seq_miss

                if c in [i[0] for i in kvc]:
                    kvv = kvc[([i[0] for i in kvc].index(c))][1]

                vertical = max(0, matrix[r-1][c] - (float(self.seq_gap)*(kvv-1) +
                        self.seq_opening))

                horizontal = max(0, matrix[r][c-1] - (float(self.seq_gap)*(kh-1) +
                        self.seq_opening))

                matrix[r][c] = sorted([diag, vertical, horizontal, 0], reverse=True)

                ## Checks gap continuation dynamically, through tuple list.

                if (vertical or horizontal) >= diag and sum([vertical,
                    horizontal, diag]) > 0:
                    if vertical >= horizontal:
                        if c in [i[0] for i in kvc]:
                            kv.append((c, kvc[([i[0] for i in kvc].index(c))][1]+1))
                        else:
                            kv.append((c, 1))

                        if horizontal == vertical:
                            kh += 1

                    else:
                        kh += 1

                else:
                    kh = 1

            kvc = kv
            print(kvc)
            kv = []

        ## Note: Implement method for users to specify prioritising matches over
        ## continuation of gap.

        return matrix


    def _pathMatrix(self, matrix):
        """Determines path through sequence scoring matricies."""

        seq1L = list(self.seq1)
        seq2L = list(self.seq2)
        outList = []

        r, c = np.unravel_index(matrix.argmax(), matrix.shape)

        while r and c != 0:
            diag = matrix[r-1][c-1]
            horizontal = matrix[r][c-1]
            vertical = matrix[r-1][c]

            g = max(diag, horizontal, vertical)

            if g == diag:
                outChar = (seq1L[r-1], seq2L[c-1])
                r -= 1
                c -= 1

            elif horizontal == vertical:
                if random.random() < 0.5:
                    outChar = (seq1L[r-1], "-")
                    r -= 1

                else:
                    outChar = ("-", seq2L[c-1])
                    c -= 1

            elif g == horizontal:
                outChar = ("-", seq2L[c-1])
                c -= 1

            else:
                outChar = (seq1L[r-1], "-")
                r -= 1

            outList.append(outChar)

        ## Condition for when r/c becomes zero needed.

        return outList


    def _seqAssembly(self, seqList):
        """Builds sequence alignment from path prediction."""

        cList = [[x[0] for x in seqList], "", [x[1] for x in seqList]]
        mList = []

        for i in range(len(cList[0])):
            if seqList[i][0] == seqList[i][1]:
                mList.append("|")
            else:
                mList.append(" ")

        cList[1] = mList

        return cList


    def _seqStatistics(self, cList):
        """Generate basic statistics on sequence alignment"""

        n = SequenceStatistics(cList)

        print(n.gc_ratio())





class SequenceStatistics:
    def __init__(self, seqList):
        self.seqList = seqList

    def gc_ratio(self):
        seq1 = (float(self.seqList[0].count("G")) +
                self.seqList[0].count("C"))/len(self.seqList[0])

        seq2 = (float(self.seqList[2].count("G")) +
                self.seqList[2].count("C"))/len(self.seqList[2])

        return seq1, seq2






x = LocalAlignment(seqA, seqB)
x.start()
