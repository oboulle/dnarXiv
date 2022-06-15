#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:56:59 2020

@author: oboulle
"""

import sys
import random
import os


def create_sequence(size: str, h_max: int) -> str:
    """
    :param size: size of the sequence to create
    :param h_max: maximum size of the homopolymer in the generated sequence
    :return: a sequence of the given size with no homopolymer larger than h_max
    """
    sequence = ""
    while len(sequence) < size:
        alphabet = ["A", "G", "C", "T"]
        letter = random.choice(alphabet)
        if sequence[-h_max:] == letter*h_max:  # if the end of the sequence is an homopolymer of this letter
            alphabet.remove(letter)
            letter = random.choice(alphabet)  # then pick another one
        sequence += letter
    return sequence


# =================== main ======================= #
if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("usage : sequence_generator.py fileName nbrSequence sizeSequence h_max")
        sys.exit(1)
    
    fileName = sys.argv[1] # file to save the sequences
    nbrSequence = int(sys.argv[2]) # number of sequences
    sizeSequence = int(sys.argv[3]) # size of the sequences
    h_max = int(sys.argv[4]) # maximum size of the homopolymers
    
    if os.path.isfile(fileName):
        print("file already exists : "+fileName)
        sys.exit(1)
    
    f = open(fileName, "w+")
    
    for i in range(nbrSequence):
        sequence = create_sequence(sizeSequence, h_max)
        f.write(">sequence_"+str(i)+"\n")
        f.write(sequence+"\n")
    f.close()
