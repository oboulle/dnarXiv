#!/usr/bin/python3

import argparse
import random
import dna_file_reader as dfr


"""
simulation of the synthesis of fragments
"""

def synthesis(input_path: str, output_path: str, nbr_synth: int, i_error: float, d_error: float, s_error: float) -> None:
    """
    create a file containing a set of duplicated fragments with errors
    """
    input_sequences = dfr.read_fasta(input_path)
    file_output = open(output_path, "w+")
    num_frag = 0
    for fragment_name, fragment in input_sequences.items():
        for i in range(nbr_synth):
            synth_line = add_errors(fragment, i_error, d_error, s_error)
            file_output.write(">frag_" + fragment_name + "_"+str(i)+ "\n")
            file_output.write(synth_line + "\n")
        num_frag += 1
    file_output.close()


def add_errors(line: str, i_error: float, d_error: float, s_error: float) -> str:
    """
    add the errors to the sequence
    :param line: line corresponding to the sequence
    :param i_error: insertion error rate
    :param d_error: deletion error rate
    :param s_error: substitution error rate
    :return: the line containing the errors
    """
    synt_line = ""
    for letter in line:
        if random.random() < s_error:  # substitution error
            alphabet = ["A", "G", "C", "T"]
            alphabet.remove(letter)
            synt_line += random.choice(alphabet)
        else:
            synt_line += letter
        if random.random() < i_error:  # insertion error
            synt_line += random.choice(["A", "G", "C", "T"])
        if random.random() < d_error:  # deletion error
            synt_line = synt_line[:-1]
    return synt_line


# =================== main =======================#
if __name__ == '__main__':
    # --------- argument part -------------#
    parser = argparse.ArgumentParser(description='generate multiple samples of each fragments with errors')
    parser.add_argument('-i', action='store', dest='input_path', required=True,
                        help='the input fasta file')
    parser.add_argument('-o', action='store', dest='output_path', required=True,
                        help='the output fasta file')
    parser.add_argument('-n', action='store', dest='nbr_synth', required=True,
                        type=int, help='number of synthesis for each fragment')
    parser.add_argument('--i_error', action='store', dest='i_error',
                        type=float, help='insertion error rate', default=0)
    parser.add_argument('--d_error', action='store', dest='d_error',
                        type=float, help='deletion error rate', default=0)
    parser.add_argument('--s_error', action='store', dest='s_error',
                        type=float, help='substitution error rate', default=0)

    print("fragment synthesis...")
    # ---------- input list ---------------#
    arg = parser.parse_args()
    synthesis(arg.input_path, arg.output_path, arg.nbr_synth, arg.i_error, arg.d_error, arg.s_error)
    print("\tcompleted !")

