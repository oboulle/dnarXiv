#!/usr/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 14:03:59 2020

@author: oboulle
"""
import random
import os
import sys
import inspect
import argparse


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir+"/synthesis_modules")
import dna_file_reader as dfr


def select_molecules(input_path: str, output_path: str, start_primer: str, stop_primer: str, n_mol: int) -> None:
    """
    randomly select molecules from the input file and print the result in the output file
    :param input_path: .fasta file containing the molecules to select
    :param output_path: .fasta file containing the selected molecules
    :param start_primer: start primer of the document to extract
    :param stop_primer: stop primer of the document to extract
    :param n_mol: number of molecule selections
    """
    molecules_dict = dfr.read_fasta(input_path)
    #filter the molecules with the corresponding start and stop primers
    corresponding_molecules = {name: seq for name, seq in molecules_dict.items()} #TODO if seq.startswith(start_primer) or seq.startswith(stop_primer)}
    pool_list = list(corresponding_molecules.items()) #list of couples (name, sequence) of molecules with the correct primers

    file_output = open(output_path, "w+")
    
    for i in range(n_mol):  
        if pool_list == []: #empty list
            print("warning molecule selection : empty pool list")
            file_output.close()
            return
        random_index = random.randint(0, len(pool_list)-1) # choice of a random index in the molecule pool
        mol_name, mol_value = pool_list[random_index] # get the associated couple

        if random.random() < .5:
            file_output.write(">" + mol_name + "\n")
            file_output.write(mol_value + "\n")
        else: # 50% chance to get the reverse complement of the molecule
            file_output.write(">rev_" + mol_name + "\n")
            file_output.write(dfr.reverse_complement(mol_value) + "\n")
        #del pool_list[random_index] #delete the couple from the pool to avoid picking it another time
    file_output.close()


# =================== main =======================#
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='selection of molecules')
    parser.add_argument('-i', action='store', dest='input_path', required=True,
                        help='.fasta file to read the molecules')
    parser.add_argument('-o', action='store', dest='output_path', required=True,
                        help='.fasta file to save the selected molecules')
    parser.add_argument('--start', action='store', dest='start_primer', required=True,
                        type=str, help='start primer of the document to extract')
    parser.add_argument('--stop', action='store', dest='stop_primer', required=True,
                        type=str, help='stop primer of the document to extract')
    parser.add_argument('-n', action='store', dest='n_mol', required=True,
                        type=int, help='number of molecules to select')

    # ---------- input list ---------------#
    arg = parser.parse_args()

    print("molecule selection...")
    select_molecules(arg.input_path, arg.output_path, arg.start_primer, arg.stop_primer, arg.n_mol)  
    print("\tcompleted !")
