#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 04 14:53:59 2020

@author: oboulle
"""
import argparse
import shutil
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0, parentdir+"/synthesis_simulation")

import utils.dna_file_reader as dfr


class FastqSequence:
    """
    class to represent a fastq sequence
    """
    def __init__(self, fastq_list):
        self.name = fastq_list[0]
        self.value = fastq_list[1]
        self.reverse_complement = dfr.reverse_complement(self.value)
        self.description = fastq_list[2]
        self.score = fastq_list[3]
        self.primers_points = {}  # 1 point is added for a primer for each kmer the sequence contains
        self.rc_primers_points = {}  # 1 point is added for a primer for each kmer the reverse complement of the sequence contains


def demultiplexing(input_path, output_dir_path, primers_dir_path, kmer_size, point_threshold):
    """
    link the fastq sequences to the corresponding primers
    for each kmer of a primer the sequence contains, there will be 1 point for this primer
    the primer with the most points will be the linked primer (if above the point threshold)
    :param input_path: .fastq file containing the sequences to demultiplex
    :param output_dir_path: directory with .fastq files for the primers containing the linked sequences
    :param primers_dir_path: directory containing the primers used in the synthesis
    :param kmer_size: size of the sub sequences of the primers used
    :param point_threshold: minimal number of kmer needed to consider a sequence linked to a primer
    """
    primers_dict = {}  # dictionary containing all the primers
    for file in os.listdir(primers_dir_path):
        primers_dict.update(dfr.read_fasta(primers_dir_path + "/" + file))

    sorted_seq_by_primers = {}  # dictionary that will contain the FastqSequences sorted by the corresponding primers
    sorted_seq_by_primers["unlinked"] = []  # add a list to save the unlinked sequences

    for fq in dfr.read_fastq(input_path):
        fastq_seq = FastqSequence(fq)
        for primer_name, primer_value in primers_dict.items():
            primer_size = len(primer_value)
            # choose the part of the sequence to compare with the primer, depending if it's a left or right primer
            if "left" in primer_name:
                # look at the beginning of the fastq_sequence and its reverse complement for left_part primers
                fastq_part = fastq_seq.value[:primer_size + 5]
                rc_fastq_part = fastq_seq.reverse_complement[:primer_size + 5]
                primer_group = primer_name.replace("left_", "")
            else:
                # look at the end of the fastq_sequence and its reverse complement for right_part primers
                fastq_part = fastq_seq.value[-primer_size - 5:]
                rc_fastq_part = fastq_seq.reverse_complement[-primer_size - 5:]
                primer_group = primer_name.replace("right_", "")
            # init a point counter for this primer in the sequence object
            if primer_group not in fastq_seq.primers_points:
                fastq_seq.primers_points[primer_group] = 0
                fastq_seq.rc_primers_points[primer_group] = 0
            # search if the kmers of the primer are in the sequence or in it's reverse complement
            for k in range(primer_size - kmer_size + 1):
                kmer = primer_value[k:k + kmer_size]
                if kmer in fastq_part:
                    fastq_seq.primers_points[primer_group] += 1
                if kmer in rc_fastq_part:
                    fastq_seq.rc_primers_points[primer_group] += 1

        # define the primer with the highest score
        primer_with_most_points = max(fastq_seq.primers_points, key=fastq_seq.primers_points.get)
        primer_points = fastq_seq.primers_points[primer_with_most_points]

        # define the primer with the highest score for the reverse complement of the sequence
        rc_primer_with_most_points = max(fastq_seq.rc_primers_points, key=fastq_seq.rc_primers_points.get)
        rc_primer_points = fastq_seq.rc_primers_points[rc_primer_with_most_points]

        # if the sequences or its reverse complement has not enough point, the sequence is not linked
        if max(primer_points, rc_primer_points) < point_threshold:
            sorted_seq_by_primers["unlinked"].append(fastq_seq)
            continue

        if primer_points > rc_primer_points:
            linked_primer = primer_with_most_points
        else:
            # the reverse complement has more points, it means that the wrong dna chain was read by the sequencer
            linked_primer = rc_primer_with_most_points
            fastq_seq.value = fastq_seq.reverse_complement

        # add the sequence in the linked primer list
        if linked_primer in sorted_seq_by_primers:
            sorted_seq_by_primers[linked_primer].append(fastq_seq)
        else:
            sorted_seq_by_primers[linked_primer] = [fastq_seq]

    try:
        os.mkdir(output_dir_path)
    except OSError:
        shutil.rmtree(output_dir_path)
        os.mkdir(output_dir_path)

    # write the fastq files containing the linked sequences for each primer
    for primer_name, fastq_sequences_list in sorted_seq_by_primers.items():
        sequences_output = open(output_dir_path + "/" + primer_name + ".fastq", "w+")
        print(primer_name, ":", len(fastq_sequences_list), "lectures associées")
        for fastq_sequence in fastq_sequences_list:
            sequences_output.write(fastq_sequence.name + "\n")
            sequences_output.write(fastq_sequence.value + "\n")
            sequences_output.write(fastq_sequence.description + "\n")
            sequences_output.write(fastq_sequence.score + "\n")
        sequences_output.close()


# =================== main =======================#
if __name__ == '__main__':
    # --------- argument part -------------#
    parser = argparse.ArgumentParser(description='sort the sequences by their primers')
    parser.add_argument('-i', action='store', dest='input_path', required=True,
                        help='the input fastq file')
    parser.add_argument('-o', action='store', dest='output_dir_path', required=True,
                        help='the output directory')
    parser.add_argument('-p', action='store', dest='primers_dir_path', required=True,
                        help='the directory of the used primers')
    parser.add_argument('--kmer_size', action='store', dest='kmer_size',
                        type=int, help='size of the subsequences of the primers to search in the fastq sequences',
                        default=10)
    parser.add_argument('--point-threshold', action='store', dest='point_threshold',
                        type=int,
                        help='threshold of the minimum number of kmer found in the fastq_sequences to link it to a primer',
                        default=10)

    # ---------- input list ---------------#
    arg = parser.parse_args()
    demultiplexing(arg.input_path, arg.output_dir_path, arg.primers_dir_path, arg.kmer_size, arg.point_threshold)
