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

if len(sys.argv) != 3:
    print("usage : result_analysis.py initial_sequences_path result_sequences_path")
    sys.exit(1)

initial_sequences_path = sys.argv[1]  # .fasta file containing the initial sequences used in the workflow
result_sequences_path = sys.argv[2]  # directory containing the resulting sequences of the workflow

initial_sequences = dfr.read_fasta(initial_sequences_path)

for seq_name, init_seq_value in initial_sequences.items():
    result_seq = dfr.read_fasta(result_sequences_path+"/"+seq_name+".fasta")
    print(seq_name, ":")
    print("init :   ", init_seq_value)
    print("result : ", result_seq["consensus_200"])
