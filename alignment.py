#!/usr/bin/python2.7

# Author: Jonathan Hull (jonathan.hull11@gmail.com)

import sys
#import threading
#import multiprocessing as mp
#import Queue
import numpy as np
import copy


seqA = "AGTAAACCGTA"
seqB = "AGTAACGTA"

class LocalAlignment:

    def __init__(self,
                seq1, 
                seq2,
                score_method="affine",
                seq_opening=5,
                seq_gap=1,
                seq_match=3, 
                seq_miss=3):

        self.seq1 = list(seq1)
        self.seq2 = list(seq2)
        self.score_method = score_method.lower()
        self.seq_match = seq_match
        self.seq_miss = seq_miss
        self.seq_gap = seq_gap
        self.seq_opening = seq_opening


    def start(self):
        """Initialise local sequence alignment algorithms"""

        self.boolMat = self.boolMatrixGen()

        if self.score_method == "affine":
            self.scoreMat = self.AffineScoreMatrix(self.boolMat)
        elif self.score_method == "linear":
            print("here")
            self.scoreMat = self.scoreMatrix(self.boolMat)
        else:
            raise NameError('score_method argument must be either '
                    '"affine, or "linear".')

        pathDict = self._pathMatrix2(self.scoreMat)
        self.base = self._seqAssembly(pathDict)


    def boolMatrixGen(self):
        """Generates boolean similarity matrix from nucleotide sequences"""

        matrix = np.zeros([len(self.seq1), len(self.seq2)], dtype=bool)

        for r in range(len(self.seq1)):
            for c in range(len(self.seq2)):
                if self.seq1[r] == self.seq2[c]:
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
                    diag = matrix[r-1][c-1] - self.seq_miss

                if c in [i[0] for i in kvc]:
                    kvv = kvc[([i[0] for i in kvc].index(c))][1]

                vertical = max(0, matrix[r-1][c] - (float(self.seq_gap)*(kvv-1) +
                        self.seq_opening))

                horizontal = max(0, matrix[r][c-1] - (float(self.seq_gap)*(kh-1) +
                        self.seq_opening))

                matrix[r][c] = sorted([diag, vertical, horizontal, 0],
                        reverse=True)[0]

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
            kv = []

        ## Note: Implement method for users to specify prioritising matches over
        ## continuation of gap.

        return matrix


    def _pathMatrix2(self, matrix):
        """Sequence alignment controller method."""

        ## Modifying list whilst iterating... (fix)
        ## Rename "outList" 

        r, c = np.unravel_index(matrix.argmax(), matrix.shape)

        indexMat = [{"r":int(r), "c":int(c), "outList":[], "sumScore":0}]

        complete_alignments = []

        ## indexMax[path] can be "both" or "vertical", pathing always goes
        ## horizontal first.

        while len(indexMat) > 0:

            indexMatc = indexMat[0]
            del(indexMat[0])

            indexMatc = self._pathMatrixMan(indexMatc)

            if 0 in [indexMatc["r"], indexMatc["c"]]:
                complete_alignments.append(indexMatc)
                continue

            elif len(indexMatc) > 1:
                indexMat.extend(indexMatc)
                continue

            else:
                indexMat.append(indexMatc)

        return complete_alignments


    def _pathMatrixMan(self, indexMat):
        """Determines paths through sequence alignment matrix."""

        r, c = indexMat["r"], indexMat["c"]

        while (0 not in [r, c]):
            
            diag = self.scoreMat[r-1][c-1]
            horizontal = self.scoreMat[r][c-1]
            vertical = self.scoreMat[r-1][c]
            indexMat["sumScore"] += self.scoreMat[r][c]

            g = [diag, horizontal, vertical]

            if g.count(max(g)) > 1:
                if horizontal == vertical > diag:
                    indexMat["sumScore"] -= self.scoreMat[r][c]
                    indexMat["r"], indexMat["c"] = r, c

                    return indexMat

                if (diag > vertical == horizontal) or g.count(0) == 3:
                    outChar = (self.seq1[r-1], self.seq2[c-1])
                    r -= 1
                    c -= 1

                elif max(g) == horizontal:
                    outChar = ("-", self.seq2[c-1])
                    c -= 1

                else:
                    outChar = (self.seq1[r-1], "-")

            elif max(g) == horizontal:
                outChar = ("-", self.seq2[c-1])
                c -= 1

            elif max(g) == vertical:
                outChar = (self.seq1[r-1], "-")
                r -= 1

            else:
                outChar = (self.seq1[r-1], self.seq2[c-1])
                r -= 1
                c -= 1

            indexMat["outList"].append(outChar)
            indexMat["r"], indexMat["c"] = int(r), int(c)

        return indexMat


    def _seqAssembly(self, pathDict):
        """Builds sequence alignment from path prediction."""

        output = []

        for adict in pathDict:
            cList = [[x[0] for x in adict["outList"]], "", [x[1] for x in adict["outList"]]]
            mList = []

            for i in range(len(cList[0])):
                if cList[0][i] == cList[2][i]:
                    mList.append("|")
                else:
                    mList.append(" ")

            cList[1] = mList

            output.append(cList)

        return output
