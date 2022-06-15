#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 09:13:17 2020

@author: oboulle
"""
import random
import inspect
import os
import sys

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import utils.dna_file_reader as dfr


def overlap_fragment(input_path, output_path, fragment_size, overlap_size):
    """
    read the input sequences and divides it into fragments with overlap
    :param input_path: path of the sequence to use (.fasta)
    :param output_path: fasta file containing the resulting fragments
    :param fragment_size: size of the fragments
    :param overlap_size: size of the overlap
    """
    input_sequences = dfr.read_fasta(input_path)
    file_output = open(output_path, "w+")
    for seq_name, seq_value in input_sequences.items():
        fragment_list = sequence_fragmentation_fixed_fragment_and_overlap(seq_value, fragment_size, overlap_size)
        i = 0
        for fragment in fragment_list:
            file_output.write(">" + seq_name + "_part" + str(i) + "\n")
            file_output.write(fragment + "\n")
            i += 1
    file_output.close()


def sequence_fragmentation_fixed_overlap(sequence, fragment_min_size, overlap_size):
    """
    divides a sequence into fragments with a fixed overlap size and variable fragments size
    :param sequence: the input sequence
    :param fragment_min_size: minimal size of the fragments
    :param overlap_size: size of the overlap
    :return: a list of the fragments
    """
    # number of fragments in the sequence and rest
    n_fragments, rest = divmod(len(sequence)-overlap_size, fragment_min_size-overlap_size)
    # the rest is distributed in the fragments
    additionnal_frag_size, rest_frag = divmod(rest, n_fragments)
    fragment_list = []
    overlap_list = []
    index = 0
    for i in range(n_fragments):
        fragment_size = fragment_min_size + additionnal_frag_size
        if i < rest_frag:
            fragment_size += 1
        fragment = sequence[index:index+fragment_size]
        fragment_list.append(fragment)
        index += fragment_size-overlap_size
        overlap = sequence[index:index+overlap_size]
        if i < n_fragments-1: overlap_list.append(overlap)
    # test if each overlap only appears once
    for overlap in overlap_list:
        if overlap_list.count(overlap) > 1:
            print("Warning, same overlap appeared multiple times ("+overlap+")")
    return fragment_list


def sequence_fragmentation_fixed_fragment(sequence, fragment_size, overlap_min_size):
    """
    divides a sequence into fragments with a fixed fragments size and variable overlap size
    :param sequence: the input sequence
    :param fragment_size: size of the fragments
    :param overlap_min_size: minimal size of the overlap
    :return: a list of the fragments
    """
    # number of fragments in the sequence and rest
    n_fragments, rest = divmod(len(sequence)-overlap_min_size, fragment_size-overlap_min_size)
    # number of letters needed for an other fragment, this number is distributed in the overlaps
    if rest == 0:
        lack = 0
        n_fragments -= 1
    else:
        lack = fragment_size-overlap_min_size-rest
    # additional_overlap_size : number of letters added in all the overlaps
    additional_overlap_size, rest_overlap = divmod(lack, n_fragments)
    fragment_list = []
    overlap_list = []
    index = 0
    print(n_fragments, rest, lack, additional_overlap_size, rest_overlap)
    for i in range(n_fragments+1):
        overlap_size = overlap_min_size + additional_overlap_size
        if i < rest_overlap:
            overlap_size += 1
        print(index,index+fragment_size,"overlap size :",overlap_size)
        fragment = sequence[index:index+fragment_size]
        fragment_list.append(fragment)
        # prepare the index for the next fragment with the corrected overlap
        index += fragment_size-overlap_size
        overlap = sequence[index:index + overlap_size]
        if i < n_fragments: overlap_list.append(overlap)
    # test if each overlap only appears once
    for overlap in overlap_list:
        if overlap_list.count(overlap) > 1:
            print("Warning, same overlap appeared multiple times ("+overlap+")")
    return fragment_list


def sequence_fragmentation_fixed_fragment_and_overlap(sequence, fragment_size, overlap_size):
    """
    divides a sequence into fragments with a fixed fragments size and overlap size, complete the last fragment if necessary
    :param sequence: the input sequence
    :param fragment_size: size of the fragments
    :param overlap_size: size of the overlap
    :return: a list of the fragments
    """
    # number of fragments in the sequence and rest
    n_fragments, rest = divmod(len(sequence)-overlap_size, fragment_size-overlap_size)
    if rest > 0:
        fragment_complement_size = fragment_size-overlap_size-rest
        for k in range(fragment_complement_size):
            sequence += random.choice(["A","C","G","T"])
        n_fragments += 1
    fragment_list = []
    overlap_list = []
    index = 0
    for i in range(n_fragments):
        fragment = sequence[index:index+fragment_size]
        fragment_list.append(fragment)
        # prepare the index for the next fragment with the corrected overlap
        index += fragment_size-overlap_size
        overlap = sequence[index:index + overlap_size]
        if i < n_fragments-1: overlap_list.append(overlap)
    # test if each overlap only appears once
    for overlap in overlap_list:
        if overlap_list.count(overlap) > 1:
            print("Warning, same overlap appeared multiple times ("+overlap+")")
    return fragment_list

    
# =================== main ======================= #
if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("usage : overlap_fragmentation.py sequences_path result_sequences_path fragment_size overlap_size")
        sys.exit(1)

    sequences_path = sys.argv[1]  # .fasta file containing the sequences used
    result_sequences_path = sys.argv[2]  # directory containing the resulting sequences
    fragment_size = int(sys.argv[3])  # size of the fragments
    overlap_size = int(sys.argv[4])  # size of the overlap

    overlap_fragment(sequences_path, result_sequences_path, fragment_size, overlap_size)
