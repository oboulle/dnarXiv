#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 09:16:59 2020

@author: oboulle
"""
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0, parentdir+"/synthesis_simulation/utils")

import dna_file_reader as dfr
import alignment_NW


def select_sequences(initial_sequences_path, result_sequences_path):
    """
    read the initial sequences and the associated results
    :param initial_sequences_path: path to the .fasta file containing the initial sequences
    :param result_sequences_path: path to the directory containing the result sequences
    :return: 2 dictionary containing the sequences
    """
    initial_sequences = dfr.read_fasta(initial_sequences_path)
    result_sequences = dict()
    size_variation = [0, 1, -1, 2, -2, 3, -3]  # consecutive variations of the size to find the closest sequence size
    for seq_name, init_seq_value in initial_sequences.items():
        result_seq_dict = dfr.read_fasta(result_sequences_path + "/" + seq_name + ".fasta")
        sequence_size = len(init_seq_value)
        result_sequences[seq_name] = None
        if result_seq_dict is None:
            continue
        for sv in size_variation:
            result_seq_name = "consensus_" + str(sequence_size + sv)
            if result_seq_name in result_seq_dict:
                result_sequences[seq_name] = result_seq_dict[result_seq_name]
                break
    return initial_sequences, result_sequences


def count_error_rates(initial_sequences, result_sequences):
    """
    count the rate for the errors in the result sequences
    :param initial_sequences: dictionary of the initial sequences
    :param result_sequences: dictionary of the result sequences
    :return:
    """
    tot_i_errors = 0
    tot_d_errors = 0
    tot_s_errors = 0
    tot_size = 0
    for seq_name, init_seq_value in initial_sequences.items():
        result_seq_value = result_sequences[seq_name]
        if result_seq_value is None: continue
        alignment_A, alignment_B = alignment_NW.alignment_nw(init_seq_value, result_seq_value)
        i_errors, d_errors, s_errors = alignment_NW.calculate_errors(alignment_A, alignment_B)
        tot_i_errors += i_errors
        tot_d_errors += d_errors
        tot_s_errors += s_errors
        tot_size += len(init_seq_value)
        print(seq_name)

    print("insertion errors :", tot_i_errors/tot_size, "; deletion errors :", tot_d_errors/tot_size, "; substitution errors :", tot_s_errors/tot_size)


if len(sys.argv) != 3:
    print("usage : result_analysis.py initial_sequences_path result_sequences_path")
    sys.exit(1)

initial_sequences_path = sys.argv[1]  # .fasta file containing the initial sequences used in the workflow
result_sequences_path = sys.argv[2]  # directory containing the resulting sequences of the workflow

initial_sequences, result_sequences = select_sequences(initial_sequences_path, result_sequences_path)
count_error_rates(initial_sequences, result_sequences)
