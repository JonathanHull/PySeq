#!/usr/bin/python2.7

# Author: Jonathan Hull (jonathan.hull11@gmail.com)

# y = threading.Thread(name="one", target=self.method)

# self.queue_one = Queue.Queue()
# self.queue_one.put(Value)
# self.queue_one.get()


import sys
import random
import threading
#import multiprocessing as mp
import Queue
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
            self.scoreMat = self.scoreMatrix(self.boolMat)
        else:
            raise NameError('score_method argument must be either '
                    '"affine, or "linear".')


        pathTup = self._pathMatrix2(self.scoreMat)

        print(pathTup)

        print("".join([x[0] for x in pathTup[0]["outList"]]))
        print("".join([x[1] for x in pathTup[0]["outList"]]))


        quit()

        pathTup = self._pathMatrix(self.scoreMat)
        base = self._seqAssembly(pathTup)

        #self._seqStatistics(base)



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

        ## use self.scoreMat intead of matrix.
        ## Only accept indexMat appendages.
        ## complete_alignments self object?

        ## Modifying list whilst iterating... (fix)
        ## Rename "outList" 

        r, c = np.unravel_index(matrix.argmax(), matrix.shape)

        indexMat = [{"r":int(r), "c":int(c), "outList":[], "path":"diagonal",
                "sumScore":0}]

        complete_alignments = []

        ## indexMax[path] can be "both" or "vertical", pathing always goes
        ## horizontal first.

        while len(indexMat) > 0:

            indexMatc = indexMat[0]
            del(indexMat[0])

            indexMatc = self._pathMatrixMan(indexMatc)

            if 0 in [indexMatc["r"], indexMatc["c"]]:
                complete_alignments.append(indexMatc)

            else:
                indexMat.extend(indexMatc)

        return complete_alignments


    def _pathMatrixMan(self, indexMat):
        """Determines paths through sequence alignment matrix."""

        ## Note: Horizontal pathing comes before vertical.
        ## Can just return the next iteration of indexMat.

        ## Note: The vert/horiz priority only makes sence if the nucleotides
        ## match, otherwise the diag should take priority.

        r, c = indexMat["r"], indexMat["c"]

        if self.scoreMat[r-1][c] == self.scoreMat[r][c-1]:
            indexMat["outList"].append((self.seq1[r-1], self.seq2[c-1]))
            indexMat["sumScore"] += self.scoreMat[r][c]
            indexMat["path"] = "both"
            return indexMat

        elif indexMat["path"] == "both":
            diagPath = indexMat
            vertPath = IndexMat
            diagPath["c"] -= 1
            vertPath["r"] -= 1
            return diagPath vertPath

        #elif indexMat["path"] == "horizontal":
            ## Need to deal with horiz/vert and pass two dicts back.
            ## Return two copies of the value with horiz/vert paths
            ## Could impletement this in the first if statement of the while
            ## loop.

            #indexMat["sumScore"] += self.scoreMat[r][c]
            #outChar = ("-", seq2L[c-1])
            #c -= 1

        #else:
        #    indexMat["sumScore"] += self.scoreMat[r][c]
        #    outChar = (self.seq1[r-1], "-")
        #    r -= 1

        while (0 not in [r, c]) and (self.scoreMat[r-1][c] !=
                self.scoreMat[r][c-1]):

            diag = self.scoreMat[r-1][c-1]
            horizontal = self.scoreMat[r][c-1]
            vertical = self.scoreMat[r-1][c]
            indexMat["sumScore"] += self.scoreMat[r][c]

            g = [diag, horizontal, vertical]

            if g.count(max(g)) > 1:
                if horizontal == vertical:
                    indexMat["sumScore"] -= self.scoreMat[r][c]
                    indexMAt["path"] == "both"
                    indexMat["r"], indexMat["c"] = r, c

                    return indexMat

                elif max(g) == horizontal:
                    outChar = ("-", self.seq2[c-1])
                    c -= 1

                else:
                    outChar = (self.seq1[r-1], "-")
                    r -= 1

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

    def _pathMatrix(self, matrix):
        """Determines path through sequence scoring matricies."""

        ## Enable multithreading
        ## Create method of storing locations where paths diverge.
        ## Create new thread where paths diverge.
        ## End state where r/c == 0, or score == 0.
        ## ^ maybe search x nucleotides after score == 0?

        outList = []

        r, c = np.unravel_index(matrix.argmax(), matrix.shape)

        while r and c != 0:
            diag = matrix[r-1][c-1]
            horizontal = matrix[r][c-1]
            vertical = matrix[r-1][c]

            g = max(diag, horizontal, vertical)

            if g == diag:
                outChar = (self.seq1[r-1], self.seq2[c-1])
                r -= 1
                c -= 1

            elif horizontal == vertical:
                if random.random() < 0.5:
                    outChar = (self.seq1[r-1], "-")
                    r -= 1

                else:
                    outChar = ("-", self.seq2[c-1])
                    c -= 1

            elif g == horizontal:
                outChar = ("-", self.seq2[c-1])
                c -= 1

            else:
                outChar = (self.seq1[r-1], "-")
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
## vertically and horizontally. Sync threads and rebuild matrix then move
# reference.
