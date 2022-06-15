import os
import sys
import consensus
import time
import inspect
from os.path import isfile, join

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0, parentdir+"/synthesis_simulation")

import utils.dna_file_reader as dfr
import utils.dna_numbering as dnbr


def compute_consensus(frag_dir):
    """
    return the list of consensus fragments from the clusters fragments directory
    """
    CONSENSUS = []
    k = 0
    
    for file in os.listdir(frag_dir):
        if isfile(join(frag_dir, file)):       
            F = []
            sequences = dfr.read_fasta(join(frag_dir, file))
            for seq_name, sequence in sequences.items():
                F.append(sequence)
     
            seqc = consensus.fragConsensus(F,FRAG_SIZE-2*len(SPACER),SPACER)
            CONSENSUS.append(seqc)
            print ("  ",k,len(seqc), seqc)
            k += 1
    return CONSENSUS


def extract_frag_numbers(CONSENSUS):
    """
    find the fragment number in each consensus fragment
    return a list of the original fragments where index = fragment number
    """
    print("Assemble fragments")
    sorted_frag_list = []
    
    for tagged_fragment in CONSENSUS:
        fragment, number = dnbr.extract_number(tagged_fragment, FRAG_SIZE, TAG_SIZE-2)
        if number >= 2*len(CONSENSUS) or number < 0:
            print("Warning, frag number out of range :",number)
        else: 
            if number >= len(sorted_frag_list):
                for i in range(number - len(sorted_frag_list)):
                    sorted_frag_list.append("_" * (FRAG_SIZE-2*len(SPACER)-TAG_SIZE))
                sorted_frag_list.append(fragment)
            else:
                sorted_frag_list[number] = fragment
        
    return sorted_frag_list


def consensus_from_2_frag_dir(path_dir_1, path_dir_2):
    """
    compute the consensus from 2 independent clusterings, 
    merge the results of the 2 consensus
    """
    CONSENSUS_1 = compute_consensus(output_dir+"/frag1")
    CONSENSUS_2 = compute_consensus(output_dir+"/frag2")
    sorted_frag_list_1 = extract_frag_numbers(CONSENSUS_1)
    sorted_frag_list_2 = extract_frag_numbers(CONSENSUS_2)
    
    seq = ""
    expected_frag_size = FRAG_SIZE-2*len(SPACER)-TAG_SIZE
    for i in range(min(len(sorted_frag_list_1), len(sorted_frag_list_2))):
        frag_1 = sorted_frag_list_1[i]
        frag_2 = sorted_frag_list_2[i]
        if frag_1 == frag_2: #same fragment
            seq += frag_1
            continue
        if frag_1.count("_") > 0: #case the frag is empty and has been filled with _ (= "_____...")
            seq += frag_2
            continue
        if frag_2.count("_") > 0:
            seq += frag_1 
            continue
        #the 2 frags are differents but not empty, the frag with the closest size to the expected size is chosen 
        if abs(len(frag_1)-expected_frag_size) > abs(len(frag_2)-expected_frag_size):
            seq += frag_2
        else:
            seq += frag_1
    return seq
    
    
if len(sys.argv) != 5:
    print("usage : reconstruct.py output_dir spacer.fasta fragment_size tag_size")
    sys.exit(1)

output_dir = sys.argv[1]
SPACER = dfr.read_fasta(sys.argv[2])["spacer"]
print ("Spacer:  ", SPACER)
FRAG_SIZE = int(sys.argv[3])
print ("Fragment size:", FRAG_SIZE)
TAG_SIZE = int(sys.argv[4])
print ("Tag size:  ", TAG_SIZE)


print ("Compute consensus fragments")
seq = consensus_from_2_frag_dir(output_dir+"/frag1", output_dir+"/frag2")


print(seq)
ff = open(output_dir+"/reconstructed_sequence.fasta","w")
ff.write (">sequence built :\n"+seq+"\n")
ff.close()
print (" - sequence size =",len(seq))


