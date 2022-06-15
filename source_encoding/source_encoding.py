#!/usr/bin/python3

import sys
import os
import inspect
import argparse
import subprocess

import hashing
import image_conversion as img_conv
import binary_conversion as bin_conv

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0, os.path.dirname(currentdir)+"/synthesis_modules")
import dna_file_reader as dfr
import dna_numbering as dnbr

"""
encode a document into a .fasta file containing the dna fragments
"""

DNA_NUMBER_LENGTH = 11 # constant in the workflow #length allocated to the dna index


def encode_document(doc_path: str, output_path: str, assembly_type: str, fragment_length: int) -> (str, str):
    """
    convert a document into a dna sequence
    :param doc_path: path to the document to encode
    :param file_extension: type of document (txt, png, etc...)
    :return: dna sequence and associated encoding data
    """
    if not os.path.isfile(doc_path):
        print("error in source encoding, input file not found at :",doc_path)
        exit(1)

    head, doc_name = os.path.split(doc_path)
    file_extension = doc_name.split(".")[1]
    
    # fragments are shorter for indexed fragmentation because the index takes some bases
    if assembly_type == "indexed":
        fragment_length -= DNA_NUMBER_LENGTH
    
    if file_extension == "png":
        sequence, encoding_data = img_conv.encode_png(doc_path, fragment_length, "RGB")
        # add a filter to the whole sequence to shuffle every bases
        sequence = apply_dna_filter(sequence, "0")
    else:
        # for binary encoding, the filter is directly applied to the binary string
        sequence, encoding_data = bin_conv.encode_file(doc_path, fragment_length), {}
        
    dfr.save_sequence_to_fasta("source", sequence, output_path) # save the source sequence
    
    encoding_data["doc_name"] = doc_name
    encoding_data["doc_type"] = file_extension

    return encoding_data


def indexed_sequence_fragmentation(sequence_path: str, total_fragment_length: int) -> list:
    """
    indexed fragmentation, add an index at the beginning of each fragments #TODO
    """
    
    _, sequence = dfr.read_single_sequence_fasta(sequence_path)
    
    
    payload_length = total_fragment_length - overhang_length - DNA_NUMBER_LENGTH
    
    # number of fragments in the sequence and rest
    n_fragments, r = divmod(len(sequence), payload_length)
    if r > 0: # if there is some rest, add an other fragment, and fill it with the start of the sequence (good/bad idea ?)
        n_fragments += 1
        sequence += sequence[:payload_length-r]
    
    indexed_fragments_list= []
    index = 0
    for num_frag in range(n_fragments):      
        dna_num = dnbr.int_to_dna_number(num_frag, DNA_NUMBER_LENGTH)
        fragment = sequence[index:index+payload_length]
        indexed_fragments_list.append(dna_num+fragment)
        # prepare the index for the next fragment
        index += payload_length
        
    return indexed_fragments_list


def apply_dna_filter(sequence: str, hash_key: str) -> list:
    """
    apply a filter to the dna sequence    
    :param sequence: dna sequence
    :param hash_key: key used to create the hash filter
    :return: the filtered sequence
    """
    base_4_to_dna_dict = {0: "A", 1: "C", 2: "G", 3: "T"}
    dna_to_base_4_dict = {"A": 0, "C": 1, "G": 2, "T": 3}
    
    filter = hashing.hash_string_to_formated_base4(hash_key, len(sequence))

    filtered_sequence = ""
    for i in range(len(sequence)):
        fragment_number = dna_to_base_4_dict[sequence[i]]     
        filter_number = int(filter[i])
        sum = (fragment_number+filter_number) % 4
         
        filtered_sequence += base_4_to_dna_dict[sum]
    return filtered_sequence


def update_metadata_file(container_path: str, doc_index: str, encoding_data: dict) -> None:
    """
    add the encoding data informations of the document to the metadata file
    """
    metadata_manager_path = os.path.dirname(currentdir)+"/workflow_commands/metadata_manager.sh"
    for key, value in encoding_data.items():
        update_command = '. '+metadata_manager_path+' && add_doc_param '+container_path+' '+doc_index+' '+key+' '+str(value)
        subprocess.call('/bin/bash -c "$ADDPARAM"', shell=True, env={'ADDPARAM': update_command})


# =================== main ======================= #
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='encode a document into dna fragments')
    parser.add_argument('-i', action='store', dest='document_path', required=True,
                        help='path to the document to encode')
    parser.add_argument('-o', action='store', dest='output_path', required=True,
                        help='path to the output fasta fragments')
    parser.add_argument('-l', action='store', dest='fragment_length', required=True,
                        type=int, help='length of the fragments')
    parser.add_argument('-t', action='store', dest='assembly_type', required=True,
                        help='assembly type of the encoded document (indexed/ordered)')
    parser.add_argument('--cont_doc', action='store', dest='container_doc_path', required=False,
                        help='path to the document in the container')

    # ---------- input list ---------------#
    arg = parser.parse_args()

    print("source encoding...")    
    
    encoding_data = encode_document(arg.document_path, arg.output_path, arg.assembly_type, arg.fragment_length)
    
    if arg.container_doc_path is not None:
        container_path, doc_index = os.path.split(os.path.normpath(arg.container_doc_path))
        update_metadata_file(container_path, doc_index, encoding_data)

    print("\tcompleted !")

