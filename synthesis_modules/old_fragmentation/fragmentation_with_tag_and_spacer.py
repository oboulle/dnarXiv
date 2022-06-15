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
import utils.dna_numbering as dnbr


def fragment(input_path, output_path, spacer_path, fragment_size, tag_size):
    """
    reads the input sequence and divides it into fragments
    :param input_path: path of the sequence to use (.fasta)
    :param output_path: fasta file containing the resulting fragments
    :param fragment_size: size of the fragments
    :param overlap_size: size of the overlap
    """
    seq_name, seq_value = dfr.read_single_sequence_fasta(input_path)
    _, spacer = dfr.read_single_sequence_fasta(spacer_path)
    
    base_fragment_list = sequence_division(seq_value, fragment_size-tag_size-2*len(spacer))
    formatted_fragment_list = fragment_formatting(base_fragment_list, spacer, tag_size)
    i = 0
    file_output = open(output_path, "w+")
    for fragment in formatted_fragment_list:
        file_output.write(">" + seq_name + "_part" + str(i) + "\n")
        file_output.write(fragment + "\n")
        i += 1
    file_output.close()


def sequence_division(sequence, base_fragment_size):
    """
    divides a sequence into fragments with a fixed fragments size, complete the last fragment if necessary
    :param sequence: the input sequence
    :param base_fragment_size: size of the fragments
    :return: a list of the fragments
    """
    if(base_fragment_size <= 0):
        print("Error : bad base fragment size :",base_fragment_size)
    # number of fragments in the sequence and rest
    n_fragments, rest = divmod(len(sequence), base_fragment_size)
    if rest > 0:
        fragment_complement_size = base_fragment_size-rest
        for k in range(fragment_complement_size):
            sequence += random.choice(["A","C","G","T"])
        n_fragments += 1
    fragment_list = []
    index = 0
    for i in range(n_fragments):
        fragment = sequence[index:index+base_fragment_size]
        fragment_list.append(fragment)
        # prepare the index for the next fragment
        index += base_fragment_size
    return fragment_list


def fragment_formatting(base_fragment_list, spacer, tag_size):
    """
    add a tag to the fragments of the list
    :param base_fragment_list: the input list of base fragments
    :param spacer: the spacer used
    :param tag_size: size of the tag
    :return: a list of the formatted fragments
    """
    nucleotides = ["A", "C", "G", "T"]
    formatted_fragment_list = []
    for i in range(len(base_fragment_list)):
        fragment = base_fragment_list[i]
        dna_number = dnbr.int_to_dna_number(i, tag_size-2)
        
        # create a nucleotide to separate the spacer from the fragment
        choices = [n for n in nucleotides if n not in (spacer[-1], fragment[0])]
        start_tag = choices[0]
        # create a nucleotide to separate the spacer from the tag
        choices = [n for n in nucleotides if n not in (dna_number[-1], spacer[0])]
        end_tag = choices[0]
        
        formatted_fragment = spacer + start_tag + fragment + dna_number + end_tag + spacer
        formatted_fragment_list.append(formatted_fragment)
        
    return formatted_fragment_list
    
    
# =================== main ======================= #
if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("usage : fragmentation.py sequences_path result_sequences_path spacer_path fragment_size tag_size")
        sys.exit(1)

    sequences_path = sys.argv[1]  # .fasta file containing the sequences used
    result_sequences_path = sys.argv[2]  # directory containing the resulting sequences
    spacer_path = sys.argv[3]  #path of the spacer file (.fasta)
    fragment_size = int(sys.argv[4])  # size of the fragments
    tag_size = int(sys.argv[5])  # size of the tag

    fragment(sequences_path, result_sequences_path, spacer_path, fragment_size, tag_size)
