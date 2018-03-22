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

seqA = "AGTAAACCGTA"
seqB = "AGTAACGTA"
seqC = "AGGTACGTATA"

class Muscle:

    def __init__(self,
            kmer_length=3,
            kmer_measure="default",
            *sequences):

        """
        Muscle Alignment

        :param kmer_length: K-mer lengths used to determine sequence similarity.
        :param kmer_measure: K-mer similarity scoring method : args - default/binary.
        :param sequences: list of sequences to be aligned.
        """

        self.kmer_length = kmer_length
        self.sequences = sequences
        self.kmer_measure = kmer_measure

        if "U" in sequences[0]:
            self.RNA = True
        else:
            self.RNA = False

    def start(self):
        """Initialises MUSCLE alignment algorithms."""

        ## Serialise k-mer dicts?

        self.kmer_base = self.kmer_dictionary()
        self.kmer_dict = self.kmer_count()

        if self.kmer_measure == "binary":
            self.kmer_binary_sim()

        else:
            self.kmer_similarity()


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
            adict["seq_{}".format(seqID+1)] = self.kmer_base

            sequence = self.sequences[seqID]

            for i in range(self.kmer_length, len(sequence)+1):
                adict["seq_{}".format(seqID+1)][sequence[i-self.kmer_length:i]]+=1

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

                    outScore += n

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


    def distance_measure(self):
        pass









x = Muscle(3,"binary",seqA,seqB, seqC)
x.start()

print(x.kmer_binary_sim())
print(x.kmer_similarity())
