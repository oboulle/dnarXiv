#!/usr/bin/python3

import os
import sys
import inspect
import subprocess
import time
import argparse


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0, os.path.dirname(currentdir)+"/synthesis_modules")
import dna_file_reader as dfr

# used to get the consensus sequence from a .fastq reads file of an ordered assembly

KMER_SIZE = 30 # arbitrary constant, must be greater than the start_sequence


def pre_filter(sequence: str) -> bool:
    """
    filter the strange sequences, that are likely sequencing errors
    True if bad sequence, False otherwise
    """
    # ridiculously long sequences
    if len(sequence) > 5000: #TODO careful, with a lot of fragments we may get very long sequences
        #print("removed too long : ",sequence, len(sequence))
        return True
    
    bases = ["A", "C", "G", "T"]

    # size h homopolymeres
    h = 15
    if any(base*h in sequence for base in bases):
        #print("removed homopolymere : ",sequence)
        return True
    
    # size h2 pairs of bases repetitions that clearly are bugs
    h2 = 10
    if any((base1+base2)*h2 in sequence for base1 in bases for base2 in bases):
        #print("removed 2-repetition : ",sequence)
        return True
      
    return False


def count_kmers_dsk(input_path: str, start_sequence: str) -> (str, dict):
    """
    use the dsk tool to count kmer
    """
    dsk_script_path = currentdir+"/dsk/build/bin/dsk"
    dsk_tempfile_path = currentdir+"/dsk/tmp/kmer_count"
    
    min_occ = 1 # threshold of occurrences to save the kmer

    kmer_count_command = dsk_script_path+' -file '+input_path+' -out '+dsk_tempfile_path+' -kmer-size '+ str(KMER_SIZE)+' -abundance-min '+str(min_occ)+' -verbose 0'
    subprocess.call('/bin/bash -c "$DSK"', shell=True, env={'DSK': kmer_count_command})
    
    # convert the result into a txt file
    dsk_ascii_script_path = currentdir+"/dsk/build/bin/dsk2ascii" # -file dsk_tempfile_path -out test.txt
    dsk_tempfile_txt_path = currentdir+"/dsk/tmp/kmer_count.txt"
    
    count_convert_command = dsk_ascii_script_path+' -file '+dsk_tempfile_path+'.h5 -out '+dsk_tempfile_txt_path+' -verbose 0'
    subprocess.call('/bin/bash -c "$CONVERT"', shell=True, env={'CONVERT': count_convert_command})
    
    kmer_occurrences_dict = {}
    
    count_file = open(dsk_tempfile_txt_path)
    line = count_file.readline()
    while line != "":
        kmer_sequence = line.split(" ")[0]
        occurrences = int(line.split(" ")[1])
        kmer_occurrences_dict[kmer_sequence] = occurrences
        kmer_occurrences_dict[dfr.reverse_complement(kmer_sequence)] = occurrences
        line = count_file.readline()
    
    potential_starters_dict = dict((k,v) for k,v in kmer_occurrences_dict.items() if k.startswith(start_sequence))
    
    if len(potential_starters_dict) == 0:
        print("kmer_consensus error : starting sequence not found in the reads")
        exit(1)
    start_kmer = max(potential_starters_dict, key=potential_starters_dict.get) # kmer with the maximum occurrences that begins with the starting sequence 
    return start_kmer, kmer_occurrences_dict


def build_sequence_from_kmer_occurrences(start_kmer: str, kmer_occurrences_dict: dict, expected_sequence_size: int) -> str:
    """
    find a list of overlapping kmers from the dict of kmer occurrences
    build the original sequence to match the expected size
    
    """
    overlap_kmers_dict = {start_kmer : True} # dict of overlapping kmers
    total_path = start_kmer # builded sequence
    
    current_kmer = start_kmer
        
    while len(total_path) < expected_sequence_size:
        
        best_next_kmer = ""
        best_next_weight = 0
        
        for next_base in ["A","C","G","T"]:
            next_kmer = current_kmer[1:]+next_base
            
            next_weight = kmer_occurrences_dict.get(next_kmer, 0)
            if next_weight == 0:
                next_weight = kmer_occurrences_dict.get(dfr.reverse_complement(next_kmer), 0)
                
            if next_weight > best_next_weight:
                best_next_kmer = next_kmer
                best_next_weight = next_weight
        
        # no best base found
        if best_next_kmer == "":
            print("kmer_consensus warning : not enough reads for the full reconstruction")
            return total_path
        
        # the chosen next base is placed at the right of the result sequence 
        total_path = total_path + best_next_kmer[-1]
        current_kmer = best_next_kmer

        # the builded sequence is starting to loop -> very bad -> need to increase the kmer size
        if overlap_kmers_dict.get(best_next_kmer, False):
            print("! loop !",best_next_kmer, best_next_weight)
            print(total_path)
            return total_path
        else:
            overlap_kmers_dict[best_next_kmer] = True
            
        #print(added_kmer, kmer_occurrences_dict.get(current_kmer,{}))
        
    return total_path


def kmer_consensus(input_path: str, output_path: str, start_sequence: str, expected_length: int) -> None:
    """
    get the consensus sequence from an ordered fragments assembly
    """
    
    #start = time.time()
    
    start_kmer, kmer_occurrences_dict = count_kmers_dsk(input_path, start_sequence)
    
    #print("\ncounting",time.time() - start,"seconds")
    #start = time.time()
    
    result_sequence = build_sequence_from_kmer_occurrences(start_kmer, kmer_occurrences_dict, expected_length) # build the resulting sequence from the dictionary of following bases
        
    #print("\nbuilding",time.time() - start,"seconds")
    
    dfr.save_sequence_to_fasta("consensus", result_sequence, output_path)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='get the original complete sequence from sequencing reads of ordered assembly')
    parser.add_argument('-i', action='store', dest='input_path', required=True,
                        help='fastq file to read the sequences')
    parser.add_argument('-o', action='store', dest='output_path', required=True,
                        help='fasta file to save the result')
    parser.add_argument('--start', action='store', dest='start_primer', required=True,
                        help='start primer at the beginning of the sequence to reconstruct')
    parser.add_argument('--stop', action='store', dest='stop_primer', required=True,
                        help='stop primer at the beginning of the reverse complement of the sequence to reconstruct')
    parser.add_argument('-e', action='store', dest='expected_length', required=True,
                        type=int, help='expected length of the output sequence')
    
    arg = parser.parse_args()
    
    print("kmer consensus...")

    kmer_consensus(arg.input_path, arg.output_path, arg.start_primer, arg.expected_length)
                
    print("\tcompleted !")

