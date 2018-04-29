#!/usr/bin/python2.7

# Author: Jonathan Hull (jonathan.hull11@gmail.com)

import sys
#import threading
#import multiprocessing as mp
#import Queue
import numpy as np
import copy

class Alignment:
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

        ## Note: holder is used so that local alignment scores do not become
        ## negative.

        matrix = np.zeros([len(boolMat)+1, len(boolMat[0])+1])
        holder = 0

        ## Initialises global scoring method.
        if self.score_method == "global":
            for c in range(len(self.seq2)+1):
                matrix[0][c] = c*self.seq_gap

            for r in range(len(self.seq1)+1):
                matrix[r][0] = r*self.seq_gap

            holder = None

        for r in range(1, len(matrix)):
            for c in range(1, len(matrix[0])):
                if boolMat[r-1][c-1] == True:
                    diag = matrix[r-1][c-1] + self.seq_match
                else:
                    diag = matrix[r-1][c-1] + self.seq_miss

                vertical = matrix[r-1][c] + self.seq_gap
                horizontal = matrix[r][c-1] + self.seq_gap

                matrix[r][c] = sorted([diag, vertical, horizontal, holder],
                        reverse=True)[0]

        return matrix

    def _seqAssembly(self, pathDict):
        """Builds sequence alignment from path prediction."""

        for adict in pathDict:
            cList = [[x[0] for x in adict["outList"]], "", [x[1] for x in adict["outList"]]]
            mList = []

            for i in range(len(cList[0])):
                if cList[0][i] == cList[2][i]:
                    mList.append("|")
                else:
                    mList.append(" ")

            cList[1] = mList

        return cList

    def _pathMatrixMan(self, matrix):
        """Sequence alignment controller method."""

        ## Modifying list whilst iterating... (fix)

        if self.score_method == "global":
            r, c = len(matrix)-1, len(matrix[0])-1

        else:
            r, c = np.unravel_index(matrix.argmax(), matrix.shape)

        indexMat = [{"r":int(r), "c":int(c), "outList":[], "sumScore":0}]
        complete_alignments = []

        while len(indexMat) > 0:
            indexMatc = indexMat[0]
            del(indexMat[0])

            indexMatc = self._pathMatrix(indexMatc)

            if 0 in [indexMatc["r"], indexMatc["c"]]:
                complete_alignments.append(indexMatc)
                continue

            elif len(indexMatc) > 1:
                indexMat.extend(indexMatc)
                continue

            else:
                indexMat.append(indexMatc)

        return complete_alignments

    def _pathMatrix(self, indexMat):
        """Determines paths through sequence alignment matrix."""

        r, c = indexMat["r"], indexMat["c"]

        while (1 not in [r, c]):
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

        ## detect end repeating nucleotide sequences which can drastically
        ## elongate alignment.

        ## Algorithm doesn't align repeat regions well; implement statistics
        ## detailing uncertainty.

        if r==1 and (self.scoreMat[r][c] >= vertical):
            outChar = (self.seq1[r-1], self.seq2[c-1])

        if c==1 and (self.scoreMat[r][c] >= horizontal):
            outChar = (self.seq1[r-1], self.seq2[c-1])

        indexMat["outList"].append(outChar)
        indexMat["sumScore"] += self.scoreMat[r][c]
        indexMat["r"] = 0

        return indexMat

    def alignment_output(self):
        return [t[::-1] for t in ["".join(x) for x in self.base]]

class LocalAlignment(Alignment):
    def __init__(self,
                seq1, 
                seq2,
                score_method="affine",
                seq_opening=-5,
                seq_gap=-1,
                seq_match=3, 
                seq_miss=-3):

        """
        Local Alignment 

        :param seq1/2: Nucleotide sequences to be aligned, no additional metadata
            should be included (i.e. fasta indentifiers).
        :param score_method: Local scoring algorithm. options are:
            - Affine
            - linear

        :param seq_opening: Nucleotide indel opening penalty. Only used with the
            Affine scoring method.
        :param seq_gap: Affine gap extension penalty / linear indel penalty.
        :param seq_match: Matching nucleotide score.
        :param seq_miss: Miss-matching nucleotide scoring penalty.

        Note: Defaults set for Affine scoring algorithm.
        Note: seq_gap = -5 better scoring option for linear method.
        """

        self.seq1 = list(seq1)
        self.seq2 = list(seq2)
        self.score_method = score_method.lower()
        self.seq_match = seq_match
        self.seq_miss = seq_miss
        self.seq_gap = seq_gap
        self.seq_opening = seq_opening

    def start(self):
        """Initialises local sequence alignment algorithms specified by user"""

        self.boolMat = self.boolMatrixGen()

        if self.score_method == "affine":
            self.scoreMat = self.AffineScoreMatrix(self.boolMat)
        elif self.score_method == "linear":
            self.scoreMat = self.scoreMatrix(self.boolMat)
        else:
            raise NameError('score_method argument must be either '
                    '"affine, or "linear".')

        pathDict = self._pathMatrixMan(self.scoreMat)
        self.score = pathDict[0]["sumScore"]
        self.base = self._seqAssembly(pathDict)
        self.alignment_list = self.alignment_output()

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
                    diag = matrix[r-1][c-1] + self.seq_miss

                if c in [i[0] for i in kvc]:
                    diag = matrix[r-1][c-1] + self.seq_miss

                if c in [i[0] for i in kvc]:
                    kvv = kvc[([i[0] for i in kvc].index(c))][1]

                vertical = max(0, matrix[r-1][c] + (float(self.seq_gap)*(kvv-1) +
                         self.seq_opening))

                horizontal = max(0, matrix[r][c-1] + (float(self.seq_gap)*(kh-1) +
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

class GlobalAlignment(Alignment):
    def __init__(self,
            seq1,
            seq2,
            seq_gap=-5,
            seq_match=3,
            seq_miss=-3):

        """
        Global Alignment 

        :param seq1/2: Nucleotide sequences to be aligned, no additional metadata
            should be included (i.e. fasta indentifiers).

        :param seq_gap: sequence indel penalty.
        :param seq_match: Matching nucleotide score.
        :param seq_miss: Miss-matching nucleotide scoring penalty.
        """

        self.seq1 = list(seq1)
        self.seq2 = list(seq2)
        self.seq_gap = seq_gap
        self.seq_match = seq_match
        self.seq_miss = seq_miss
        self.score_method = "global"

    def start(self):
        """Initialise global sequence alignment algorithms"""

        self.boolMat = self.boolMatrixGen()
        self.scoreMat = self.scoreMatrix(self.boolMat)
        pathDict = self._pathMatrixMan(self.scoreMat)
        self.score = pathDict[0]["sumScore"]
        self.base = self._seqAssembly(pathDict)
        self.alignment_list = self.alignment_output()
