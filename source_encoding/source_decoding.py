#!/usr/bin/python3

import sys
import os
import subprocess
import inspect
import argparse

import hashing
import image_conversion as img_conv
import binary_conversion as bin_conv

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0, os.path.dirname(currentdir)+"/synthesis_modules")
import dna_file_reader as dfr
import dna_numbering as dnbr

"""
decode the resulting dna fragments of the workflow to recreate the original document
"""

DNA_NUMBER_SIZE = 11 # constant in the workflow

def indexed_fragments_to_sequence(frag_path: str, frag_length: int, n_frag: int) -> str:
    """
    remove the filter to the fragments defined by a hash of the fragment number
    :param frag_path: path to a list of the fragments
    :param frag_length: length of a fragment
    :param n_frag: number of fragments in the original document
    :return: the original sequence
    """
    reconstitued_fragment_dict = {}
    
    fragment_list = list(dfr.read_fasta(frag_path).values())
    
    # get the correct index and remove the filter for each fragment
    for fragment in fragment_list:
        
        while len(fragment) < frag_length: # extend fragment if too short
            fragment += "A"
        if len(fragment) > frag_length:
            fragment = fragment[:frag_length] # cut fragment if too long
        
        # extract the index sequence from the fragment
        dna_index = fragment[:DNA_NUMBER_SIZE]
        # get the associated number
        num_frag = dnbr.dna_number_to_int(dna_index)
        if num_frag == -1:
            #print("unrecognized number :",dna_index)
            continue
        
        if num_frag in reconstitued_fragment_dict:
            print("fragment of number",num_frag,"already read")
        else:
            reconstitued_fragment_dict[num_frag] = fragment[DNA_NUMBER_SIZE:]
    
    # assembly of the fragments
    original_sequence = ""

    for i in range(n_frag):
        if i in reconstitued_fragment_dict:
            original_sequence += reconstitued_fragment_dict[i]
        else:
            #missing fragment, fill with _ instead
            original_sequence += "_"*(frag_length-DNA_NUMBER_SIZE)
                
    return original_sequence


def remove_dna_filter(sequence: str, hash_key: str) -> str:
    """
    remove the applied filter on a sequence
    the filter is a SHA256 hash in base 4 added to the sequence
    """
    if not sequence:
        return ""
    
    base_4_to_dna_dict = {0: "A", 1: "C", 2: "G", 3: "T"}
    dna_to_base_4_dict = {"A": 0, "C": 1, "G": 2, "T": 3}
    
    filter = hashing.hash_string_to_formated_base4(hash_key, len(sequence))

    unfiltered_sequence = ""
    for i in range(len(sequence)):
        sequence_number = dna_to_base_4_dict[sequence[i]]     
        filter_number = int(filter[i])
        sub = (sequence_number-filter_number) % 4
         
        unfiltered_sequence += base_4_to_dna_dict[sub]
    return unfiltered_sequence


def decode_document(sequence: str, output_path: str, encoding_data_dict: dict) -> None:
    """
    decode the sequence and save the decoded document
    """
    if encoding_data_dict["doc_type"] == "png":
        # remove the filter applied in the encoding
        original_sequence = remove_dna_filter(sequence, "0")
        img_conv.decode_png(original_sequence, output_path, encoding_data_dict)
    else:
        # the filter has already been removed in the binary decoding
        bin_conv.decode_file(sequence, output_path)
        
        
def get_encoding_data(encoding_data) -> dict:
    """
    get the data used to decode the document from the xml metadata file of the container
    """
    encoding_data_list = " ".join(encoding_data).split(";")
    
    encoding_data_dict = {}
    for couple in encoding_data_list:
        if len(couple) > 0:
            data_name, data_value = couple.split("=")
            encoding_data_dict[data_name] = data_value
    return encoding_data_dict


# =================== main ======================= #
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='decode a sequence')
    parser.add_argument('-i', action='store', dest='sequence_path', required=True,
                        help='path to the consensus sequence fasta file')
    parser.add_argument('-o', action='store', dest='output_path', required=True,
                        help='path to the output document file')
    parser.add_argument('-t', action='store', dest='assembly_type', required=True,
                        help='assembly type of the encoded document (indexed/ordered)')
    parser.add_argument('-l', action='store', dest='frag_length', required=False,
                        type=int, help='length of the dna fragments, required for indexed assembly')
    parser.add_argument('--data', action='store', dest='encoding_data', required=True,
                        nargs='+', help='encoding data for the document to decode')

    # ---------- input list ---------------#
    arg = parser.parse_args()
    
    print("source decoding...")
    
    encoding_data_dict = get_encoding_data(arg.encoding_data)
    
    if arg.assembly_type == "indexed":
        # using indexed fragments
        sequence = indexed_fragments_to_sequence(arg.sequence_path, arg.frag_length, int(encoding_data_dict["fragment_number"])) 
        # using ordered fragments already reassembled in a single sequence
    else:
        _, sequence = dfr.read_single_sequence_fasta(arg.sequence_path)
        
    sequence_wo_primers = sequence[20:-20] # remove the primers of the result

    decode_document(sequence_wo_primers, arg.output_path, encoding_data_dict)

    print("\tcompleted !")

