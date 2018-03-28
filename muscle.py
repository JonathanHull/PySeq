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

import numpy
import math
import copy
import alignment

from SeqStatistics import SeqStatistics

seqA = "AGTAAACCGTA"
seqB = "AGTAACGTA"
seqC = "AGGTACGTATA"

class Muscle:

    def __init__(self,
            kmer_length=3,
            kmer_measure="default",
            distance_measure="kimura",
            *sequences):

        """
        Muscle Alignment

        :param kmer_length: K-mer lengths used to determine sequence similarity.
        :param distance_measure: Sequence distance scoring method. Options are:
            - Kimura
            - Kmer

        :param kmer_measure: K-mer similarity scoring method. Options are:
            - default
            - binary

        :param sequences: list of sequences to be aligned.
        """

        self.kmer_length = kmer_length
        self.kmer_measure = kmer_measure
        self.distance_measure = distance_measure
        self.sequences = sequences

        if "U" in sequences[0]:
            self.RNA = True
        else:
            self.RNA = False

    def start(self):
        """Initialises MUSCLE alignment algorithms."""

        ## Serialise k-mer dicts?

        if self.distance_measure == "kmer":

            self.kmer_base = self.kmer_dictionary()
            self.kmer_dict = self.kmer_count()

            if self.kmer_measure == "binary":
                self.kmer_dict = self.kmer_binary_sim()

            else:
                self.kmer_dict = self.kmer_similarity()

            self.sequence_distance = self.distance_kmer()


        elif self.distance_measure == "kimura":
            self.sequence_distance = self.distance_kimura()


        else:
            raise NameError("distance_measure must be \"kimura\" or \"kmer\".")

        ## kmura distance method call

    def kmer_dictionary(self):
        """Generates dictionary of all possible k-mers given k-mer length."""

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

    def kmer_count(self):
        """Determines sequence k-mer occurrence."""

        adict = dict()

        for seqID in range(len(self.sequences)):
            adict["seq_{}".format(seqID+1)] = copy.deepcopy(self.kmer_base)

            sequence = self.sequences[seqID]

            for i in range(len(sequence)-self.kmer_length+1):
                adict["seq_{}".format(seqID+1)][sequence[i:i+self.kmer_length]]+=1

        return adict

    def kmer_similarity(self):

        ## F = sum(min[Nx(T), Ny(T)]) / [min(Lx,Ly) - k + 1]
        ## Nx,Ny == Occurance of specific kmer.
        ## Lx,Ly == Sequence lenghts.
        ## k == k-mer length.
        ## self.kmer_length

        kmers = self.kmer_dict["seq_1"].keys()
        similarity_scores = dict()

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

                similarity_scores[(prime_sequence+1, i+1)] = float(outScore)/d

        return similarity_scores


    def kmer_binary_sim(self):

        ## F = sum(deltaxy(T)) / [min(Lx,Ly) - k +1]
        ## deltaxy(y) = 1 if kmer present in both sequences.
        ## Lx,Ly == Sequence lenghts.
        ## k == k-mer length.
        ## self.kmer_length

        kmers = self.kmer_dict["seq_1"].keys()
        similarity_scores = dict()

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

                similarity_scores[(prime_sequence+1, i+1)] = float(outScore)/d

        return similarity_scores

    def distance_kimura(self):

        ## dKimura = -log(e)(1-D - D^2/5)

        distance_measure = dict()

        for first in range(len(self.sequences)):
            for second in range(first+1, len(self.sequences)):

                GA = alignment.GlobalAlignment(self.sequences[first],
                        self.sequences[second])

                GA.start()
                GA = GA.base
                D = SeqStatistics(GA).fractionalIdentity()

                distance_measure[(first+1, second+1)] = -math.log(1-(D-D**2/5))

        return distance_measure

    def distance_kmer(self):

        distance_measure = dict()

        for i in self.kmer_dict:
            distance_measure[i] = 1 - self.kmer_dict[i]

        return distance_measure


x = Muscle(3,"default","kmer",seqA,seqB,seqC)
x.start()
