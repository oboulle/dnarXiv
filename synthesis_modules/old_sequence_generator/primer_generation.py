#!/usr/bin/python3

import sys
import os
import inspect
import random
import argparse
import subprocess

import dna_file_reader as dfr


"""
choose a couple of primers for each document in the container
"""


def generate_primers(container_path: str, primers_list_path: str) -> None:
    """
    choose the primers for the given container from an existing primer list
    each document of the container receive 2 associated primers
    :return: nothing, the primers are saved in the container metadata file
    """
    path, dirs, files = next(os.walk(container_path))
    nbr_documents = len(dirs)
    primer_list = [] # generated primers
    
    fragments_path = os.path.join(container_path, "container_fragments.fasta")
    
    container_fragments = dfr.read_fasta(fragments_path)
    
    index_table_dict = {} # each key is a kmer, the associated value is the list of container fragments containing this kmer
    kmer_size = 10
    # fill the index table
    for frag_name, sequence in container_fragments.items():
        for j in range(0, len(sequence)-kmer_size+1):
            kmer = sequence[j:j+kmer_size]
            if kmer in index_table_dict:
                if frag_name not in index_table_dict[kmer]:
                    index_table_dict[kmer].append(frag_name)
            else:
                index_table_dict[kmer] = [frag_name]
    
    f = open(primers_list_path, "r")
    primers_candidates = f.readlines()
    f.close()
    
    # searching 2*nbr_documents compatible primers inside the list of potential primers
    # a compatible primer must have an hybridization rate lower than the threshold with any fragments of the container or with an other selected primer
    hybrid_threshold = 18
    while len(primer_list) < 2*nbr_documents:
        if len(primers_candidates) == 0:
            print("Error : the list of potential primers has been emptied without finding enough ("+2*nbr_documents+") primers for the container")
            exit(1)
        # select a candidate primer from the list
        candidate = primers_candidates.pop(0).replace("\n","")
        is_close_to_any_frag = False
        # test if it has common kmers with any fragments of the container or with an already selected primer
        for j in range(0, len(candidate)-kmer_size+1):
            if is_close_to_any_frag:
                break
            candidate_kmer = candidate[j:j+kmer_size]
            if candidate_kmer in index_table_dict:
                for frag_name in index_table_dict[candidate_kmer]:
                    if hybrid_rate(candidate, container_fragments[frag_name]) > hybrid_threshold:
                        # candidate kmer is too close to a fragment
                        print("too close :", candidate, "\n"+ container_fragments[frag_name],"with hybrid rate of",str(hybrid_rate(candidate, container_fragments[frag_name])))
                        is_close_to_any_frag = True
                        break
            for selected_primer in primer_list:
                if hybrid_rate(candidate, selected_primer) > hybrid_threshold:
                        # candidate kmer is too close to a selected primer
                        print("too close :", candidate, "\n", selected_primer,"with hybrid rate of",str(hybrid_rate(candidate, selected_primer)))
                        is_close_to_any_frag = True
                        break
        if not is_close_to_any_frag:
            # candidate can be selected
            primer_list.append(candidate)
    
    # save the selected primers
    for i in range(nbr_documents):
        start_primer = primer_list[2*i]
        stop_primer = primer_list[2*i+1]
        save_primer(container_path, dirs[i], start_primer, stop_primer)
        
        
def hybrid_rate(short_sequence: str, sequence: str) -> int:
    """
    give the maximum hybridization rate for a short sequence into a bigger one
    """
    max_score = 0
    for i in range(len(sequence)-len(short_sequence)+1):
        local_score = 0
        reverse_local_score = 0 # score of hybridization of the reverse complement
        for j in range(len(short_sequence)):
            if short_sequence[j] == sequence[i+j]:
               local_score += 1
            if dfr.reverse_complement(short_sequence)[j] == sequence[i+j]:
               reverse_local_score += 1
        max_score = max(max_score, local_score, reverse_local_score)
    return max_score


def save_primer(container_path: str, doc_index: str, start_primer: str, stop_primer: str) -> None:
    """
    write the primers to the container metadata file
    """
    currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # get directory of this script
    metadata_manager_path = os.path.dirname(currentdir)+"/workflow_commands/metadata_manager.sh" # link to the metadata manager script
            
    update_command = '. '+metadata_manager_path+' && add_doc_param '+container_path+' '+doc_index+' start_primer '+start_primer
    subprocess.call('/bin/bash -c "$ADDPARAM"', shell=True, env={'ADDPARAM': update_command})
    
    update_command = '. '+metadata_manager_path+' && add_doc_param '+container_path+' '+doc_index+' stop_primer '+stop_primer
    subprocess.call('/bin/bash -c "$ADDPARAM"', shell=True, env={'ADDPARAM': update_command})


def generate_random_primer_list(output_path: str, list_size: int, primer_size: int) -> None:
    file_output = open(output_path, "w+")    
    for i in range(list_size):
        primer = ""
        alphabet = ["A", "C", "G", "T"]
        for i in range(primer_size):
            primer += random.choice(alphabet)
        file_output.write(primer + "\n")
    file_output.close()


# =================== main =======================#
if __name__ == '__main__':
        
    parser = argparse.ArgumentParser(description='encode a document into dna fragments')
    parser.add_argument('-c', action='store', dest='container_path', required=True,
                        help='path to the container')

    # ---------- input list ---------------#
    arg = parser.parse_args()
    
    print("primer generation...")

    primer_list_path = os.path.realpath(__file__).replace("primer_generation.py", "primer_list.txt")#TODO path of the primer list
    generate_primers(arg.container_path, primer_list_path)  

    print("\tcompleted !")

