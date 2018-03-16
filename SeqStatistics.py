class SequenceStatistics:
    def __init__(self, seqList):
        self.seqList = seqList

    def gc_ratio(self):
        """Calculates the GC ratio of the passed nucleotide sequences."""

        seq1 = (float(self.seqList[0].count("G")) +
                self.seqList[0].count("C"))/len(self.seqList[0])

        seq2 = (float(self.seqList[2].count("G")) +
                self.seqList[2].count("C"))/len(self.seqList[2])

        return seq1, seq2

    def sequences(self):
        """Returns passed nucleoitde sequences as list."""
        return ["".join(self.seqList[0]), "".join(self.seqList[2])]

    def similarity(self):
        """Returns sequence similarity ratios."""
        return float(self.seqList[1].count("|"))/len(self.seqList[1])

    def score(self):
        """Returns sum of alignment score."""
        pass
