#!/usr/bin/python2.7

## Author: Jonathan Hull (jonathan.hull11@gmail.com)
## MUSCLE: MUltiple Sequence Comparison by Log-Expectation.

## Three Stages:
## 1 - The draft progressive

## 1.1 - Similarity measure
##     - The similarity of each pari of sequences is computed, either using
##       k-mer counting or by construcing a global alignment of the pair
##       determinig the fraactal identity.

## 2 - Improved progressive
## --- Kimura distance used to reestimate binary tree used to create draft
## --- alignment.

## 3 - Refinement stages

## Introduce method of producing heatmap for similarity matrix.

import numpy
import math
import copy
import alignment

import numpy as np
from SeqStatistics import SeqStatistics

seqA = "AGTAAACCGTA"
seqB = "AGTAACGTA"
seqC = "AGGTACGTATA"

class MuscleAlignemnt:
    def __init__(self,
            kmer_length=3,
            distance_measure="kimura",
            *sequences):

        """
        Muscle Alignment

        :param kmer_length: K-mer lengths used to determine sequence similarity.
        :param distance_measure: Sequence distance scoring method. Options are:
            - Kimura
            - Kmer
            - kmer_binary

        :param sequences: list of sequences to be aligned.
        """

        self.kmer_length = kmer_length
        self.distance_measure = distance_measure
        self.sequences = sequences

        if "U" in sequences[0]:
            self.RNA = True
        else:
            self.RNA = False


    def start(self):
        """Initialises MUSCLE alignment algorithms."""

        ## Serialise k-mer dicts?

        self.distance_matrix = self.distance_init()

        if self.distance_measure == "kimura":
            self.distance_kimura()

        elif self.distance_measure in ["kmer", "kmer_binary"]:

            self.kmer_base = self._kmerDictGen()
            self.kmer_dict = self._kmerCount()

            if self.distance_measure == "kmer_binary":
                self.kmer_binary_sim()

            else:
                self.kmer_similarity()

            self.initial_distance = self.distance_kmer()

        else:
            raise NameError("distance_measure parameter must be \"kimura\", "
                    "\"kmer\", or \"kmer_binary\".")

        self._distance_mirror()


    def _kmerDictGen(self):
        """
        _kmerDictGen

        Generates dictonary of all posible kmers of given length
        (self.kmer_length) for DNA/RNA sequences.

        i.e. Counts kmer occurance in parsed sequences.
        """

        if self.RNA == True:
            nucList = ["A","U","G","C"]
        else:
            nucList = ["A","T","G","C"]

        outlist = nucList
        kmer_count = 4**self.kmer_length
        outDict = dict()

        while len(outlist) < kmer_count:

            holder_list = []

            for i in outlist:
                for e in range(4):
                    holder_list.append(i+nucList[e])

            outlist = holder_list

        for i in outlist:
            outDict[i] = 0

        return outDict

    def _kmerCount(self):
        """
        _kmerCount

        Utilises output dictionary from _kmerDictGen to count occurance of kmer
        in parsed nucleotide sequences.
        """

        adict = dict()

        for seqID in range(len(self.sequences)):
            adict["seq_{}".format(seqID+1)] = copy.deepcopy(self.kmer_base)

            sequence = self.sequences[seqID]

            for i in range(len(sequence)-self.kmer_length+1):
                adict["seq_{}".format(seqID+1)][sequence[i:i+self.kmer_length]]+=1

        return adict


    def distance_init(self):
        """Creates basic matrix framework"""
        return np.zeros([len(self.sequences), len(self.sequences)])

    def _distance_mirror(self):
        """Mirrors distance_matrix values"""
        
        for i in range(len(self.distance_matrix)):
            for e in range(len(self.distance_matrix)):
                self.distance_matrix[e][i] = self.distance_matrix[i][e]


    def kmer_similarity(self):

        ## F = sum(min[Nx(T), Ny(T)]) / [min(Lx,Ly) - k + 1]
        ## Nx,Ny == Occurance of specific kmer.
        ## Lx,Ly == Sequence lenghts.
        ## k == k-mer length.
        ## self.kmer_length

        kmers = self.kmer_dict["seq_1"].keys()

        for prime_sequence in range(len(self.kmer_dict)):
            for i in range(prime_sequence+1, len(self.sequences)):

                outScore = 0
                    
                d = min(
                    len(self.sequences[prime_sequence]),
                    len(self.sequences[i])
                    ) - self.kmer_length + 1

                for kmer in kmers:
                    n = min(
                        self.kmer_dict["seq_{}".format(prime_sequence+1)][kmer],
                        self.kmer_dict["seq_{}".format(i+1)][kmer]
                        )

                    outScore += n

                self.distance_matrix[prime_sequence][i] = float(outScore)/d


    def kmer_binary_sim(self):

        ## F = sum(deltaxy(T)) / [min(Lx,Ly) - k +1]
        ## deltaxy(y) = 1 if kmer present in both sequences.
        ## Lx,Ly == Sequence lenghts.
        ## k == k-mer length.
        ## self.kmer_length

        kmers = self.kmer_dict["seq_1"].keys()

        for prime_sequence in range(len(self.kmer_dict)):
            for i in range(prime_sequence+1, len(self.sequences)):

                outScore = 0

                d = min(
                    len(self.sequences[prime_sequence]),
                    len(self.sequences[i])
                    ) - self.kmer_length + 1

                for kmer in kmers:
                    if min(self.kmer_dict["seq_{}".format(prime_sequence+1)][kmer],
                        self.kmer_dict["seq_{}".format(i+1)][kmer]) >= 1:

                        n = 1

                    else:
                        n=0

                    outScore += n

                self.distance_matrix[prime_sequence][i] = float(outScore)/d


    def distance_kimura(self):

        ## dKimura = -log(e)(1-D - D^2/5)

        for first in range(len(self.sequences)):
            for second in range(first+1, len(self.sequences)):

                GA = alignment.GlobalAlignment(self.sequences[first],
                        self.sequences[second])

                GA.start()
                GA = GA.base
                D = SeqStatistics(GA).fractionalIdentity()

                self.distance_matrix[first][second] = -math.log(1-(D-D**2/5))



    def distance_kmer(self):
        """Generates distance matrix utilising dkmer equation."""

        dkmer_mat = np.zeros([len(self.sequences), len(self.sequences)])

        for r in range(len(self.distance_matrix)):
            for c in range(len(self.distance_matrix[0])):
                if self.distance_matrix[r][c] == 0:
                    dkmer_mat[r][c] = 0

                else:
                    dkmer_mat[r][c] = 1 - self.distance_matrix[r][c]

        return dkmer_mat


    def __nearest_neighbour(self, distance_matrix):
        """Finds next set of sequence neighbours in distance matrix."""

        ## Find smallest value
        ## r,c,val

        sm = [0,0,1]
        holder = 1

        for row in range(len(distance_matrix)):
            for col in range(holder, len(distance_matrix[row])):

                if distance_matrix[row][col] < sm[2]:
                    sm = [row, col, distance_matrix[row][col]]

            holder+=1

        return sm

    def _alphabetic_identification(self):
        """Generates dictionary, keys are numeric cluster identifier."""

        static_list = [chr(x) for x in range(97,123)]
        sequence_ids = list(static_list)
        outlist =  list(static_list)

        while len(sequence_ids) < len(self.sequences):

            holder_list = []

            for i in outlist:
                for e in range(len(static_list)):
                    holder_list.append(i+static_list[e])

            outlist = holder_list
            sequence_ids.extend(holder_list)

        return sequence_ids[:len(self.sequences)]


    def binary_tree(self):
        """Unweighted Pair Group Method with Arithmetic Mean."""

        ## len(self.distance_matrix)-1 == number of iterations.

        test = self.distance_matrix

        for distance_iteration in range(len(self.distance_matrix)-2):
            self.distance_matrix_update(test)

            quit()


    def distance_matrix_update(self, distance_matrix):
        """Generates next iteration of the UPGMA clustering method"""

        ## column/row +1 == seq_id

        new_matrix = np.zeros([len(distance_matrix)-1, len(distance_matrix)-1])
        cluster_pair = self.__nearest_neighbour(distance_matrix)
        clustering_dict = dict()
        clustering_list = list()

        ## Use list to track clustering.
        ## Should alphabetise clustering

        ## Check whether not apart of group first
        ## then check if in tuple.
        ## [x for x in clustering_list if type(x) == tuple]

        ## dict[D1] = First clustering
        ## dict[D1] = [cluster_pair[r], cluster_pair[c], 
        ## num_clustered_seqs-1]

        #clustering_dict["D{}".format(len(clustering_dict)+1)] =

        #[cluster_pair[r], cluster_pair[c], 






















x = MuscleAlignemnt(3,"kmer",seqA,seqB,seqC)
x.start()
x.binary_tree()
