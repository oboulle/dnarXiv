#!/usr/bin/python3

import argparse
import os
import sys
import inspect
import subprocess
import random


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0, os.path.dirname(currentdir)+"/source_encoding")
import dna_numbering as dnbr

import dna_file_reader as dfr
import sequence_control as sc
import alignment_global as ag

OVERHANG_LENGTH = 0 # CONSTANT for length of overhang # also fixed in binary_conversion.encode_file() # also for expected length in dna_read
PRIMER_LENGTH = 20 # CONSTANT for length of primers

"""
design of fragments assembly with overhangs and primers :
START : start primer (20 bases)
STOP : stop primer (20 bases)
Ox : overhang (4 bases)
_ : payload

START______Oa     _________Ob  ...     _________Oj     ______STOP    (forward fragments)
START______     Oa_________    ...   Oi_________     Oj______STOP    (reverse fragments)

fragments are assembled by couple forward/reverse, then linked in correct order with the overhangs

the start primer is placed at the beginning of first fragment
the stop primer is placed at the beginning of the reverse complement of the last fragment
"""



def sequence_fragmentation(sequence_path: str, total_fragment_length: int) -> list:
    """
    split the sequence into fragments
    makes room for the 2 primers and the overhangs
    returns the list of forward fragments
    """
    
    _, sequence = dfr.read_single_sequence_fasta(sequence_path)
    
    payload_length = total_fragment_length - OVERHANG_LENGTH

    # number of fragments in the sequence and rest
    n_fragments, r = divmod(len(sequence) + 2*PRIMER_LENGTH, payload_length)
    if r > 0: # there shouldn't be rest, as last fragment is filled in binary conversion
        print("ordered fragmentation, no round number of fragments : rest",str(r))
        exit(1)
        #n_fragments += 1
        #filling_sequence = dfr.complement(sequence[:payload_length-r]) # the filling should respects the dna constraints, but not be a part of the sequence
        #sequence += filling_sequence
    
    fragment_list = []
    
    # first fragment is shorter to have place for the start primer later
    index = payload_length-PRIMER_LENGTH
    first_fragment = sequence[:index]
    
    first_fragment_wo_ban_word = remove_ban_words_z_encoding(first_fragment)
    fragment_list.append(first_fragment_wo_ban_word)
    
    for i in range(1, n_fragments-1):
        fragment = sequence[index:index+payload_length]

        # remove banned words from the fragment, the offset of the 5 bases window of the z_method used in encoding depends on the index
        fragment_wo_ban_word = remove_ban_words_z_encoding(fragment, -index % 5)
        fragment_list.append(fragment_wo_ban_word)
        
        # prepare the index for the next fragment
        index += payload_length
    
    # last fragment have the stop primer later
    last_fragment = sequence[index:]
    
    last_fragment_wo_ban_word = remove_ban_words_z_encoding(last_fragment, -index % 5)
    fragment_list.append(last_fragment_wo_ban_word)
    
    return fragment_list


def remove_ban_words_z_encoding(sequence: str, z_method_offset=0) -> str:
    """
    remove banned words from a sequence encoded with the binary_conversion.binary_to_dna_z() method
    the z method ensure a GC% in [40;60] on a window of 5 bases (&z&&z), but we the constraint is on larger windows
    changing a base to avoid a banned word should not affect the whole sequence too much
    the changed base must be a z base, and be A<=>G; T<=>C, so the decoding of the updated sequence will still return the same result
    
    if the original sequence has been fragmented, the 5 bases windows of the z method can be offset, z_method_offset is used to correct this
    """

    # change a z base but keep the binary meaning for the decoding
    dna_change_z = {"A": "G", "T": "C", "G": "A", "C": "T"}

    # check for banned words
    indexs_ban_words = []
    for banned_word in sc.get_forbidden_sequences():
        index = sequence.find(banned_word)
        while index != -1:
            indexs_ban_words.append([index, index+len(banned_word)])
            index = sequence.find(banned_word, index+1)
    
    if indexs_ban_words == []:
        return sequence
    
    
    # indexes of z bases in the sequence (2nd and 5th base every 5 base window)
    z_indexes = [ z_method_offset+ 5*k+i for k in range(-1, len(sequence)//5) for i in [1,4] ]

    sequence_wo_bans = sequence

    for banned_words_indexes in indexs_ban_words:
        # indexes of z bases that can be changed to remove a banned word
        potential_changes_index = [ k for k in z_indexes if k >= banned_words_indexes[0] and k < banned_words_indexes[1] ]

        for change_index in potential_changes_index:
            sequence_wo_bans_list = list(sequence_wo_bans) # turn into list because python strings are immutable
            
            # change a z_base without changing the meaning of encoded data
            sequence_wo_bans_list[change_index] = dna_change_z[sequence_wo_bans[change_index]]
            
            index_before_ban_word = max(0, banned_words_indexes[0]-3) # test the changed part from 3 bases before to avoid creation of homopolymeres
            index_after_ban_word = min(len(sequence)-1, banned_words_indexes[1]+3) # test to 3 bases after
            
            changed_sub_sequence = ''.join(sequence_wo_bans_list)[index_before_ban_word:index_after_ban_word]

            # test the sequence only for the banned words and homopolymeres, ignore the GC%
            if sc.sequence_check(changed_sub_sequence, window_size=60, min_GC=0, max_GC=100, verbose=False):   
                # keep the change that removes this ban word if it doesn't create homopolymeres or other ban words            
                sequence_wo_bans = ''.join(sequence_wo_bans_list)
                print("removed ban word", sequence[banned_words_indexes[0]:banned_words_indexes[1]],"->",changed_sub_sequence)
                break
            else:
                print("failed to remove ban word", sequence[banned_words_indexes[0]:banned_words_indexes[1]],"-X>",changed_sub_sequence)
                print("trying again...")
    
    # last check for very very odd cases (maybe overlapping ban words)
    for banned_word in sc.get_forbidden_sequences():
        if banned_word in sequence_wo_bans:
            print("binary conversion: unable to remove forbidden sequences", sequence_wo_bans)
            exit(1)
    
    return sequence_wo_bans

    
def select_compatible_overhangs(fragment_list: list, overhangs_path: str) -> list:
    """
    assign overhangs between each couple of fragments
    return the dict of forward and reverse fragments with their overhangs
    the fragments are placed in the dict with the following order:
    first; first_rev; 2nd; 2nd_rev; ...; last; last_rev
    """
    
    n_fragment_couple = len(fragment_list)-1

    # used when no overhangs in the assembly (only for simulations)
    if OVERHANG_LENGTH == 0:
        return ["" for k in range(n_fragment_couple)]
    
    overhang_list = list(dfr.read_fasta(overhangs_path).values())
    
    if len(overhang_list) < n_fragment_couple:
        print("error fragment design : not enough overhangs ("+str(len(overhang_list))+"), needed "+str(n_fragment_couple))
        exit(1)

    # build the compatibility matrix between fragments couple and overhangs
    # a couple is compatible if the concatenation (frag1 + overhang + frag2) respects the dna constraints
    compatibility_mat = []
    for i in range(n_fragment_couple):
        compatibility_mat.append([])
        for overhang_sequence in overhang_list:
            assembly_test = fragment_list[i] + overhang_sequence + fragment_list[i+1] # test the assembly of a fragment couple joined by the overhang
            assembly_compatibility = sc.sequence_check(assembly_test, window_size=60, min_GC=40, max_GC=60) # True if the assembly respects the dna constraints
            compatibility_mat[i].append(assembly_compatibility)
        
    
    # possible optimization : algorithme Hongrois
    def rec_get_selected_overhang(row_index=0, selected_columns=[]) -> list:
    ###
    # get the list of selected compatible overhangs in order 
    # use the compatibility_matrix, where rows are a fragment couple and columns an overhang
    # compatibility_mat[x][y] = True if the fragment couple of index x is compatible with the overhang of index y
    # return a list of overhang indexes, or [] if no solutions
    # warning, complexity is high
    ###
        # exit and end the recursion when enough overhang are selected
        if row_index == len(compatibility_mat):
            return selected_columns
    
        for i, overhang_compatibility in enumerate(compatibility_mat[row_index]):
            # try overhangs not already selected and compatible with the fragments couple
            if not i in selected_columns and overhang_compatibility:
                # go deeper in the recursion with this overhang selected
                rec_selected_columns = rec_get_selected_overhang(row_index+1, selected_columns +[i])
                # if the recursion reached a working solution, return it to previous layers until exiting the recursive function
                if rec_selected_columns != []:
                    return rec_selected_columns
                # if rec_selected_columns is empty, that means using this overhang for this couple is a dead end, so the loop continue
        
        # all possible overhangs have been tested for this row and no solution found, this is a dead end, 
        # return to the previous layer to try different overhangs 
        return []
    
    selected_overhangs_index_list = rec_get_selected_overhang()
    if selected_overhangs_index_list == []:
        # there isn't any working combination of couple fragment / overhang
        print("error in fragment design : no compatible overhangs found for the fragments")
        exit(1)
    
    selected_overhangs_list = [overhang_list[k] for k in selected_overhangs_index_list]
    return selected_overhangs_list


def add_overhangs(fragment_list: list, selected_overhangs_list: list) -> dict:
    """
    add the overhangs to the fragments list
    return a dict of forward and reverse fragments with the overhangs
    """
    fragments_with_overhangs_dict = {}
    
    for i in range(len(fragment_list)):
        
        # overhangs are placed at the end of forward fragments
        if i == len(fragment_list)-1: # no overhang at the end of the last fragment
            forward_total_fragment = fragment_list[i]
        else:
            forward_total_fragment = fragment_list[i]+selected_overhangs_list[i]
        fragments_with_overhangs_dict[str(i)+"_Fw"] = forward_total_fragment
        
        # reverse complement of overhangs are placed at the end of reverse complement fragments
        if i == 0: # no overhang at the end of the reverse complement fragment of the first fragment
            reverse_total_fragment = dfr.reverse_complement(fragment_list[i])
        else:
            reverse_total_fragment = dfr.reverse_complement(selected_overhangs_list[i-1]+fragment_list[i])
        fragments_with_overhangs_dict[str(i)+"_Rv"] = reverse_total_fragment
        
    return fragments_with_overhangs_dict


def generate_primers(fragments_with_overhangs_dict: dict) -> (str, str):
    """
    use primer_generator script to generate 2 primers that doesn't hybridize with the sequence or each others
    """
    
    #___generate a lot of primers with good properties___#
    
    primers_gen_script = currentdir+"/primer_generator/bin/primers_gen"
    random_sequence_path = currentdir+"/primer_generator/data/random_seq_10000.fasta" # some very long random sequences to create the primers from
    primers_param_path = currentdir+"/primer_generator/data/param_dnarxiv.txt" # file containing the parameters for the primer_generator scripts
    output_primers_path = currentdir+"/primer_generator/temp/primers.fasta" # generated primers
    
    # generate a file containing a lot of primers full filling the general conditions set in param_dnarxiv.txt
    primers_gen_command = primers_gen_script+' '+random_sequence_path+' '+primers_param_path+' '+output_primers_path
    subprocess.call('/bin/bash -c "$PRIMERGEN"', shell=True, env={'PRIMERGEN': primers_gen_command})
    
    
    #___select primers without hybridation with the encoded data___#
    
    total_encoded_sequence_path = currentdir+"/primer_generator/temp/total_encoded_sequence.fasta" # sequence representing the assembly of all encoded data with the overhangs
    total_sequence = ""
    fragments_with_overhangs_list = list(fragments_with_overhangs_dict.values())

    # create a file containing the concatenation of fragments with overhangs
    for i in range(0, len(fragments_with_overhangs_list), 2): # add only the forward fragments 
        total_sequence += fragments_with_overhangs_list[i]
    dfr.save_sequence_to_fasta("total_sequence", total_sequence, total_encoded_sequence_path)
    
    primers_filter_script = currentdir+"/primer_generator/bin/primers_chk"
    output_compatible_primers_path = currentdir+"/primer_generator/temp/compatible_primers.fasta" # filtered primers that doesn't hybridize with the assembly
    output_hyb = currentdir+"/primer_generator/temp/hyb" # list of hybridized primers (unused)
    output_dic = currentdir+"/primer_generator/temp/dic" # dict of primers (unused)
    output_pos = currentdir+"/primer_generator/temp/pos" # positions of primers (unused)
    
    # filter the primers that hybridize with the whole sequence of encoded data assembled with the overhangs
    filter_primers_command = primers_filter_script+' '+total_encoded_sequence_path+' '+output_primers_path+' '+primers_param_path \
                                    +' '+output_compatible_primers_path+' '+output_hyb+' '+output_dic+' '+output_pos
    subprocess.call('/bin/bash -c "$PRIMERFILTER"', shell=True, env={'PRIMERFILTER': filter_primers_command})
    
    
    #___select 2 primers that doesn't hybridize with each other___#
    
    primers_list = list(dfr.read_fasta(output_compatible_primers_path).values())
    random.shuffle(primers_list)
    
    for i in range(len(primers_list)-1):
        first_primer = primers_list[i]
        # the addition of the first primer and first forward fragment must respect dna constraints
        concatened_primer_with_first_frag = first_primer + fragments_with_overhangs_list[0]
        if not sc.sequence_check(concatened_primer_with_first_frag, window_size=60, min_GC=40, max_GC=60):
            continue
        
        for second_primer in primers_list[i+1:]:
            # the addition of the 2nd primer and last reverse complement fragment must respect dna constraints
            concatened_primer_with_last_frag = second_primer + fragments_with_overhangs_list[-1]
            if not sc.sequence_check(concatened_primer_with_last_frag, window_size=60, min_GC=40, max_GC=60):
                continue
            
            # count number of consecutive hybridized bases between the 2 primers = consecutive common bases with the complement
            max_hybridisation = ag.max_common_bases(first_primer, dfr.complement(second_primer))

            if max_hybridisation <= 4:
                return first_primer, second_primer
        
    # not a single couple of primers has been found
    print("error in fragment design : no correct primers couple found")
    exit(1)
    
    
def add_primers(fragments_with_overhangs_dict: dict, start_primer: str, stop_primer: str) -> dict:
    """
    add the primers to the first and last fragments (forward and reverse)
    """
    
    first_fragment_key = list(fragments_with_overhangs_dict)[0] 
    first_fragment_reverse_key = list(fragments_with_overhangs_dict)[1]
    
    # start primer before first fragment, rev comp of start primer after rev comp of first fragment
    fragments_with_overhangs_dict[first_fragment_key] = start_primer + fragments_with_overhangs_dict[first_fragment_key]
    fragments_with_overhangs_dict[first_fragment_reverse_key] = fragments_with_overhangs_dict[first_fragment_reverse_key] + dfr.reverse_complement(start_primer)

    last_fragment_key = list(fragments_with_overhangs_dict)[-2]
    last_fragment_reverse_key = list(fragments_with_overhangs_dict)[-1]
    
    # rev comp of stop primer after last fragment, stop primer before rev comp of last fragment
    fragments_with_overhangs_dict[last_fragment_key] = fragments_with_overhangs_dict[last_fragment_key] + dfr.reverse_complement(stop_primer)
    fragments_with_overhangs_dict[last_fragment_reverse_key] = stop_primer + fragments_with_overhangs_dict[last_fragment_reverse_key]

    return fragments_with_overhangs_dict


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
    
    parser = argparse.ArgumentParser(description='divide a sequence into fragments, choose a couple of primers and a set of overhang between the fragments')
    parser.add_argument('-i', action='store', dest='sequence_path', required=True,
                        help='path to the sequence to fragment (.fasta)')
    parser.add_argument('-o', action='store', dest='output_path', required=True,
                        help='path to the output fasta fragments')
    parser.add_argument('-l', action='store', dest='total_fragment_length', required=True,
                        type=int, help='total length of a fragment')
    parser.add_argument('--cont_doc', action='store', dest='container_doc_path', required=False,
                        help='path to the document in the container')


    # ---------- input list ---------------#
    arg = parser.parse_args()

    print("fragment design...")    
    
    # split the sequence
    fragment_list = sequence_fragmentation(arg.sequence_path, arg.total_fragment_length)
    
    # select some overhangs
    selected_overhangs_list = select_compatible_overhangs(fragment_list, currentdir+"/overhangs_10_reals.fasta")
    
    # add the overhangs to the fragments
    overhang_fragments_dict = add_overhangs(fragment_list, selected_overhangs_list)
    
    # generate 2 primers
    start_primer, stop_primer = generate_primers(overhang_fragments_dict)
        
    # add the primers to the fragments
    primer_overhang_fragments_dict = add_primers(overhang_fragments_dict, start_primer, stop_primer)
    
    # save the final fragments
    dfr.save_dict_to_fasta(primer_overhang_fragments_dict, arg.output_path)

    # optional : update the metadata of the container
    if arg.container_doc_path is not None:
        encoding_data_dict = {"fragment_number":len(fragment_list), "start_primer":start_primer, "stop_primer":stop_primer}
        container_path, doc_index = os.path.split(os.path.normpath(arg.container_doc_path))
        update_metadata_file(container_path, doc_index, encoding_data_dict)

    print("\tcompleted !")

