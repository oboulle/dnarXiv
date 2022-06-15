#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 09:26:22 2020

@author: oboulle
"""
import os
import sys


def read_fasta(fasta_file_path: str) -> dict:
    """
    :param fasta_file_path: path to the .fasta file
    :return: a dictionary containing the name and the content of the sequences in the file
    """
    sequence_dict = {}
    if not os.path.isfile(fasta_file_path):
        return sequence_dict
    
    fasta_file = open(fasta_file_path)
    line = fasta_file.readline()
    while line != "":
        if line.startswith(">"):
            sequence_name = line[1:].replace("\n", "")
            sequence = ""
            line = fasta_file.readline()
            while not line.startswith(">") and line != "":
                sequence += line.replace("\n", "")
                line = fasta_file.readline()
            sequence_dict[sequence_name] = sequence
        else:
            print("fasta format error :",line)
            exit(1)
    fasta_file.close()
    return sequence_dict


def read_single_sequence_fasta(fasta_file_path: str) -> (str, str):
    """
    :param fasta_file_path: path to the .fasta file
    :return: the name and the value of the sequence in the file
    """
    if not os.path.isfile(fasta_file_path):
        return None, None
    
    fasta_file = open(fasta_file_path)
    line = fasta_file.readline()
    if line.startswith(">"):
        sequence_name = line[1:].replace("\n", "")
        sequence = ""
        line = fasta_file.readline()
        while not line.startswith(">") and line != "":
            sequence += line.replace("\n", "")
            line = fasta_file.readline()
    else:
        print("error, could not read the sequence from the file :",fasta_file_path)
        sequence_name, sequence = None, None
    fasta_file.close()
    return sequence_name, sequence
    
    
def read_fastq(fastq_file_path: str) -> list:
    """
    :param fastq_file_path: path to the .fastq file
    :return: a list containing the name, content, description and score of the sequences in the file
    """
    sequence_list = []
    if not os.path.isfile(fastq_file_path):
        return sequence_list
    
    fastq_file = open(fastq_file_path)
    line = fastq_file.readline()
    while line != "":
        if line.startswith("@"):
            seq_name = line.replace("\n", "")
            line = fastq_file.readline()
            seq_value = line.replace("\n", "")
            line = fastq_file.readline()
            seq_description = line.replace("\n", "")
            line = fastq_file.readline()
            seq_score = line.replace("\n", "")
            line = fastq_file.readline()
            sequence_list.append([seq_name, seq_value, seq_description, seq_score])
        else:
            print("fastq format error :",line)
            exit(1)
    fastq_file.close()
    return sequence_list


def complement(sequence: str) -> str:
    """
    return the complement of a dna sequence (AGTT -> TCAA)
    :param sequence: the input sequence
    :return: the complement of the input sequence
    """
    tab = str.maketrans("ACTG", "TGAC")

    return sequence.translate(tab)


def reverse_complement(sequence: str) -> str:
    """
    return the reverse complement of a dna sequence (AGTT -> AACT)
    :param sequence: the input sequence
    :return: the reverse complement of the input sequence
    """
    tab = str.maketrans("ACTG", "TGAC")

    return sequence.translate(tab)[::-1]


def fasta_to_fastq(fasta_file_path: str, output_path: str) -> None:
    """
    convert a fasta file into a fastq file
    :param fasta_file_path: path of the .fasta file
    :param output_path: path of the .fastq output
    """
    sequences = read_fasta(fasta_file_path)
    file_output = open(output_path, "w+")
    for name, value in sequences.items():
        file_output.write("@"+name+"\n")
        file_output.write(value+"\n")
        file_output.write("+\n")
        file_output.write("---\n")
    file_output.close()
    

def fastq_to_fasta(fastq_file_path: str, output_path: str) -> None:
    """
    convert a fastq file into a fasta file
    :param fastq_file_path: path of the .fastq file
    :param output_path: path of the .fasta output
    """
    sequences = read_fastq(fastq_file_path)
    file_output = open(output_path, "w+")
    for seq_name, seq_value, seq_description, seq_score in sequences:
        file_output.write(">"+seq_name[1:]+"\n")
        file_output.write(seq_value+"\n")
    file_output.close()


def save_sequence_to_fasta(sequence_name: str, sequence: str, output_path: str) -> None:
    """
    save a single sequence in a file with the .fasta format
    """
    file_output = open(output_path, "w")
    file_output.write(">"+sequence_name+"\n")
    file_output.write(sequence+"\n")
    file_output.close()
    

def save_dict_to_fasta(sequences_dict: dict, output_path: str) -> None:
    """
    save all sequences from a dict in a file with the .fasta format
    """
    if sequences_dict is None:
        return
    
    file_output = open(output_path, "w+")
    for name, sequence in sequences_dict.items():
        file_output.write(">" + name + "\n")
        file_output.write(sequence + "\n")
    file_output.close()


# =================== main ======================= #
if __name__ == '__main__':
    input_path = sys.argv[1]  # file to read the sequences
    output_path = sys.argv[2]  # file to save the sequences
    fastq_to_fasta(input_path, output_path)
