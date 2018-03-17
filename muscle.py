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

class Muscle:

    def __init__(self,
            kmer_length=3,
            *sequences):

        """
        Muscle Alignment

        :param kmer_length: K-mer lengths used to determine sequence similarity.
        :param sequences: list of sequences to be aligned.
        """

        self.kmer_length = kmer_length
        self.sequences = sequences

        if "U" in sequences[0]:
            self.RNA = True
        else:
            self.RNA = False

    def start(self):
        """Initialises MUSCLE alignment algorithms."""

        ## Serialise k-mer dicts?

        self.kmer_base = self.kmer_dictionary()
        self.kmer_dict = self.kmer_similarity()


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


    def kmer_similarity(self):
        """Determines sequence k-mer occurrence."""

        adict = dict()

        for seqID in range(len(self.sequences)):
            adict["seq_{}".format(seqID+1)] = self.kmer_base

            sequence = self.sequences[seqID]

            for i in range(self.kmer_length, len(sequence)+1):
                adict["seq_{}".format(seqID+1)][sequence[i-self.kmer_length:i]]+=1

        return adict



x = Muscle(3,seqA,seqB)
x.start()
