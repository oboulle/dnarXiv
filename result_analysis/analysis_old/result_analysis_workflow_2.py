#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 09:16:59 2020

@author: oboulle
"""

import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0, parentdir+"/synthesis_simulation")

import utils.dna_file_reader as dfr
import utils.alignment_NW as aNW


def count_error_rates(init_seq_value, result_seq_value):
    """
    count the rate for the errors in the result sequence
    :param initial_sequence: value of the initial sequences
    :param result_sequence: value of the result sequences
    """
    alignment_A, alignment_B = aNW.alignment_nw(init_seq_value, result_seq_value)
    print(alignment_A)
    print(alignment_B)
    i_errors, d_errors, s_errors = aNW.calculate_errors(alignment_A, alignment_B)
    seq_size = len(init_seq_value)
    print("insertion errors :", i_errors/seq_size, "("+str(i_errors)+")")
    print("deletion errors :", d_errors/seq_size, "("+str(d_errors)+")")
    print("substitution errors :", s_errors/seq_size, "("+str(s_errors)+")")


if len(sys.argv) != 3:
    print("usage : result_analysis.py initial_sequences_path result_sequences_path")
    sys.exit(1)

initial_sequence_path = sys.argv[1]  # .fasta file containing the initial sequence used in the workflow 2
result_sequence_path = sys.argv[2]+"/reconstructed_sequence.fasta"  # directory containing the resulting sequence of the workflow 2

init_seq_value = next(iter(dfr.read_fasta(initial_sequence_path).values()))
result_seq_value = next(iter(dfr.read_fasta(result_sequence_path).values()))
count_error_rates(init_seq_value, result_seq_value)
