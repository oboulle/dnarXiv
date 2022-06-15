#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 09:26:22 2020

@author: oboulle
"""
import os
import sys

""""
This script converts decimal numbers into small dna sequences representing the number


First, we convert all 16 pairs of nucleotides into hexadecimal digits :
AA = 0
AC = 1
AG = 2
AT = 3
CA = 4
...
TG = e
TT = f

So we can convert our decimal number into hexadecimal and then replace each hexa digit with the associated pair of nucleotide

Problem : writing 0 in dna using all the allocated size would make AAAAAAAAAA...
while we don't want long homopolymeres or repetitive kmers

We will use a nucleotide to separate each pair of digits
0 = AAXAAXAAXAAXAA with X a nucleotide != A

In general, we will have dna numbers like N1N2 X1 N3N4 X2 N5N6 ..., 
with X1 != N2 and N3 (so 2 possible choices for X1)
and X2 != N4,N5 and X1 (so 1 possible choice for X2)

For example : 79 in decimal will be written 4f in hexa, so CA TT in dna sequence
If we allocate 8 nucleotides for the dna number, we can write :
AA  G  CA  C  TT  -> AAGCACTT
 0      4      f

This method avoids long homopolymeres and repetitions of identical kmers
"""

base_16_to_dna_dict = {"0": "AA", "1": "AC", "2": "AG", "3": "AT",
                       "4": "CA", "5": "CC", "6": "CG", "7": "CT",
                       "8": "GA", "9": "GC", "a": "GG", "b": "GT",
                       "c": "TA", "d": "TC", "e": "TG", "f": "TT"}

def int_to_dna_number(n_base_10, number_size):
    """
    convert a decimal number to the dna number 
    :param n_base_10: the number in base 10 to convert
    :param number_size: the number of nucleotides used to write the number in dna, must be 2+3*n like
    :return: the number written in nucleotides
    """
    nucleotides = ["A", "C", "G", "T"]
    n_base_16 = hex(n_base_10)[2::]
    
    dna_number_len = 3*len(n_base_16) - 1
    
    # too large number
    if dna_number_len > number_size:
        print("dna_numbering error : not enough space ("+str(number_size)+" nucleotides) to write the number",str(n_base_10))   
        exit(1)
        
    # complete the number with 0 at the beginning
    while dna_number_len < number_size:
        n_base_16 = "0"+n_base_16
        dna_number_len += 3
    
    dna_number = ""
    for digit in n_base_16:
        dna_digit = base_16_to_dna_dict[digit]
        
        if dna_number != "": 
            last_nucleotide = dna_number[-1]
        else: last_nucleotide = ""
        next_nucleotide = dna_digit[0]
        #print(dna_digit, last_nucleotide, next_nucleotide)
        # find a filler that is different from the 2 surrounding nucleotides
        filler = nucleotides.pop(0)
        while filler == last_nucleotide or filler == next_nucleotide:
            nucleotides.append(filler)
            filler = nucleotides.pop(0)
            
        nucleotides.append(filler)
        dna_number += filler + dna_digit
        
    dna_number = dna_number[1:] # remove the first filler at the beginning of the number
    return dna_number
    
    
def dna_number_to_int(dna_number):
    """
    convert a dna number into the associated decimal number
    :param dna_number: the number written in nucleotides
    :return: the number in base 10, or -1 if the dna_number was not recognized
    """
    
    base_16_number = ""
    for i in range(0, len(dna_number),3):
        dna_digit = dna_number[i:i+2]
        if dna_digit in base_16_to_dna_dict.values():
            for key, value in base_16_to_dna_dict.items():
                if value == dna_digit:
                    base_16_number += key
                    break
        else:
            return -1
        
    if base_16_number == "":
        return -1
    base_10_number = int(base_16_number, 16)
    
    #one last test, if we convert the int number into dna again, it should be strictly equal to dna_number
    if int_to_dna_number(base_10_number, len(dna_number)) != dna_number:
        #if not then the number is not recognized
        return -1
    
    return base_10_number


# =================== main ======================= #
if __name__ == '__main__':
    i=4
    #print(i, base_10_to_dna(i,5))
    for i in range(1000):
        print(i, base_10_to_dna(i,11))
