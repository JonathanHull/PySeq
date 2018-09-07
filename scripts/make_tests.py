#!/usr/bin/env python3.6

## Script used to generate data for unit testing.

import sys
import os
import pickle

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from pyseq.alignment import *

output_dir = "../tests/test_data"

def align(seq1,
          seq2, 
          alignment_type):
    """
    Align

    Generates data to run unit tests on

    Parameters:-
    :seq1/2:            Nucleotide sequence to be aligned
    :alignment_type:    Alignment algorithm. options are:
        - global
        - local
        - affine

    Note: Affine alignment is a local varient.
    Note: Make sure that the output is correct.
    """

    if alignment_type.lower() == "global":
        x = GlobalAlignment(seq1, seq2)
    elif alignment_type.lower() == "local":
        x = LocalAlignment(seq1, seq2, "linear")
    if alignment_type.lower() == "affine":
        x = LocalAlignment(seq1, seq2)

    x.start()
    return x

def save_alignment(alignment_object, output_file):
    """
    Save alignment

    Pickles the sequence alignment object.

    Parameters:-
    :alignment_object:  Alignment object (after alignment).
    :output_file:       Path to pickled file containing alignment result.
    """
    with open(output_file, 'wb') as output_file:
        pickle.dump(alignment_object, output_file)


def load_alignment(pickled_alignment):
    """
    Load alignment

    Loads pickled sequence alignment data.

    Parameters:-
    :pickled_alignment:  File path to pickled object.
    """
    with open(pickled_alignment, 'rb') as input_file:
        alignment_object = pickle.load(input_file)

    return alignment_object

test_alignments = [("AAACCCTT","AACCTTGGAA", "local"),
                   ("AAGGCATAC", "AGGGCATTCA", "local"),
                   ("ACGTAACCCCCAAAATAG", "AGCAAATCCCAAA", "local"),
                  ]

count = 1
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

for x,y,z in test_alignments:
    print("Done: {}".format(count))
    alignment_object = align(x,y,z)
    save_alignment(alignment_object, output_dir+("/{}_".format(z))+str(count))
    count += 1
