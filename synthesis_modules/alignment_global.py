#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:12:23 2020

@author: oboulle
"""

import sys
import numpy as np
import dna_file_reader as dfr


def alignment_nw(sequence_a: str, sequence_b: str) -> (str, str):
    """
    Needleman–Wunsch algorithm : find an alignment for the 2 sequences
    :param sequence_a:
    :param sequence_b:
    :return: an alignment
    """
    penalty = -1

    n_A = len(sequence_a)
    n_B = len(sequence_b)

    mat = np.zeros((n_A + 1, n_B + 1))

    for i in range(n_A + 1):
        mat[i, 0] = penalty * i
    for j in range(n_B + 1):
        mat[0, j] = penalty * j

    for i in range(n_A):
        for j in range(n_B):
            similarity = 1 if sequence_a[i] == sequence_b[j] else -1
            match = mat[i, j] + similarity
            delete = mat[i, j + 1] + penalty
            insert = mat[i + 1, j] + penalty
            mat[i + 1, j + 1] = max(match, delete, insert)

    alignment_A = ""
    alignment_B = ""
    i = n_A
    j = n_B
    while i > 0 or j > 0:
        score = mat[i, j]
        score_diag = mat[i - 1, j - 1]
        score_left = mat[i - 1, j]
        similarity = 1 if sequence_a[i - 1] == sequence_b[j - 1] else -1
        if i > 0 and j > 0 and score == score_diag + similarity:
            alignment_A = sequence_a[i - 1] + alignment_A
            alignment_B = sequence_b[j - 1] + alignment_B
            i -= 1
            j -= 1
        elif i > 0 and score == score_left + penalty:
            alignment_A = sequence_a[i - 1] + alignment_A
            alignment_B = "_" + alignment_B
            i -= 1
        else:
            alignment_A = "_" + alignment_A
            alignment_B = sequence_b[j - 1] + alignment_B
            j -= 1
    return alignment_A, alignment_B


def count_errors(base_sequence: str, copied_sequence: str) -> (int, int, int):
    """
    count the number of errors in the copied sequence
    :param base_sequence:
    :param copied_sequence:
    :return: number of insertion, deletion and substitution errors
    """
    i_errors = 0
    d_errors = 0
    s_errors = 0
    for i in range(len(base_sequence)):
        if base_sequence[i] == "_":
            i_errors += 1
        elif copied_sequence[i] == "_":
            d_errors += 1
        elif base_sequence[i] != copied_sequence[i]:
            s_errors += 1
    return i_errors, d_errors, s_errors

    
def max_common_bases(sequence_A: str, sequence_B: str):
    """
    return the maximum number of consecutive bases in the global alignment of 2 sequences
    """
    alignment_A, alignment_B = alignment_nw(sequence_A, sequence_B)
    max_common = 0
    current_common = 0
    for i in range(len(alignment_A)):
        if alignment_A[i] == alignment_B[i]:
            current_common += 1
            max_common = max(max_common, current_common)
        else:
            current_common = 0

    return max_common
   
    
# =================== main ======================= #
if __name__ == '__main__':

    if len(sys.argv) != 3:
        print("usage : alignment_global.py base_sequence other_sequence")
        sys.exit(1)

    sequence_A = sys.argv[1]
    sequence_B = sys.argv[2]

    max_common_bases(sequence_A, dfr.complement(sequence_B))
    

