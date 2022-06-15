import sys
import math


import dna_file_reader as dfr

"""
methods to control if a sequence verify some dna constraints
"""


forbidden_sequences = ["GGTCTC", "GAGACC"] # constant, sequences that must be absent from any sequence, the reverse complement must also be in the list


def get_homopolymere_nbr(sequence: str, max_h: int) -> int:
    """
    count the number of homopolymeres larger than h_max in the sequence
    """
    h_nbr = 0 #number of homopolymere larger than h_max found
    row_size = 0 #size of the current row of consecutive nucleotides
    last_nucleotide = "" #previous nucleotide in the sequence
    for nucleotide in sequence:
        if nucleotide == last_nucleotide:
            row_size += 1
        else:
            if row_size > max_h:
                h_nbr += 1
            row_size = 1
        last_nucleotide = nucleotide
    if row_size > max_h:
        h_nbr += 1
        
    return h_nbr


def get_GC_percent(sequence: str, window_size: int) -> (float, float):
    """
    returns the minimum and maximum GC percentage for all the windows of the sequence
    """
    
    def check_GC_Window(window):
        A_content = window.count('A')
        C_content = window.count('C')
        T_content = window.count('T')
        G_content = window.count('G')
        
        GC_percent = round(100*(G_content + C_content) / (A_content + C_content + T_content + G_content), 1)
        return GC_percent
    
    if len(sequence) <= window_size:
        GC_percent = check_GC_Window(sequence)
        return GC_percent, GC_percent
        
    max_GC_percent = 0
    min_GC_percent = 100
    
    for index in range(len(sequence)-window_size):
        window = sequence[index:index+window_size]
        GC_percent = check_GC_Window(window)
        max_GC_percent = max(GC_percent, max_GC_percent)
        min_GC_percent = min(GC_percent, min_GC_percent)
        
    return min_GC_percent, max_GC_percent


def get_loop_nbr(sequence: str, loop_size: int, window_size: int) -> int:
    """
    count the number of potential loops (reverse complement of a sub sequence in a local window) in the sequence
    """
    loop_nbr = 0 #number of loop
    if len(sequence) < 2*loop_size: #sequence too short for any loop
        return 0
    
    for index in range(len(sequence)-loop_size):
        sub_sequence = sequence[index:index+loop_size]
        rev_compl_sub_sequence = dfr.reverse_complement(sub_sequence)
        if rev_compl_sub_sequence in sequence[index+loop_size:index+loop_size+window_size]:
            loop_nbr += 1

    return loop_nbr


def get_forbidden_sequences_nbr(sequence: str) -> bool:
    """
    check if some forbidden sequences are in the sequence
    """
    seq_nbr = 0
    for seq in forbidden_sequences:
        if seq in sequence:
            seq_nbr += 1
    return seq_nbr
    

def sequence_check(sequence: str, window_size=60, h_size=3, min_GC=45, max_GC=55, verbose=False) -> bool:
    """
    test if a the conditions for a correct sequence are met, return True if all 4 constraints are valid
    """
    h_nbr = get_homopolymere_nbr(sequence, h_size)
    if verbose: print("number of homopolymere larger than",h_size,":",h_nbr)
    
    min_GC_percent, max_GC_percent = get_GC_percent(sequence, window_size)
    if verbose: print("GC percentage :",min_GC_percent,"% to",max_GC_percent,"%")
    
    loop_nbr = get_loop_nbr(sequence, 11, window_size)
    if verbose: print("number of potential loop :",loop_nbr)
    
    forbidden_sequences_nbr = get_forbidden_sequences_nbr(sequence)
    if verbose: print("number of forbidden_sequences :",forbidden_sequences_nbr)
    
    if h_nbr == 0 and min_GC_percent >= min_GC and max_GC_percent <= max_GC and loop_nbr == 0 and forbidden_sequences_nbr == 0:
        if verbose: print("sequence is correct")
        return True
    else:
        if verbose: print("sequence is not correct")
        return False
    
    
def get_forbidden_sequences() -> list:
    """
    just return the list of forbidden sequences
    """
    return forbidden_sequences
    
    
def hash_until_correct(sequence: str, start_key=0) -> (int, str):
    """
    use a hash on the sequence until it meets the conditions, return the hash key and the hashed sequence
    """
    hash_key = start_key
    hashed_sequence = source_encoding.apply_filter(sequence, str(hash_key))
    while not sequence_check(hashed_sequence):
        hash_key += 1
        hashed_sequence = source_encoding.apply_filter(sequence, str(hash_key))
    
    return hash_key, hashed_sequence
    
    
def find_hash_keys(sequence: str) -> None:
    """
    find the keys to hash the sequence to pass the conditions
    the sequence is divided in sub sequences and each one is hashed until it passes the conditions
    
    """
    tot_hash = ""
    hash_keys = []
    sub_seq_size = 150 #higher -> less keys, but more time consuming 
    sub_seq_nbr = int(math.ceil(len(sequence)/sub_seq_size))
    for i in range(sub_seq_nbr):
        #display progress bar
        k=int(20*len(hash_keys)/sub_seq_nbr)
        sys.stdout.write('\r')
        sys.stdout.write("[%-20s] %d%%" % ('='*k, 5*k))
        sys.stdout.flush()
        
        sub_seq = sequence[sub_seq_size*i:sub_seq_size*(i+1)]
        
        hash_key, hashed_sub_seq = hash_until_correct(sub_seq)
        
        #if the total hashed sequence do not pass the conditions, the sub_seq is re hashed to find an other key
        while not sequence_check(tot_hash+hashed_sub_seq):
            hash_key, hashed_sub_seq = hash_until_correct(sub_seq, hash_key+1)
        
        
        hash_keys.append(hash_key)
        tot_hash += hashed_sub_seq
    
    print("\n"+tot_hash)
    print(hash_keys)
    sequence_check(tot_hash, 60, True) #will be valid


def decode_sequence(sequence: str) -> None:
    hash_keys = [74, 9257, 4300, 1690, 6164, 276, 3452, 4847, 333, 541, 835, 6263, 31, 29842, 904, 856, 3397, 3044, 360, 690, 273, 459, 5217, 13906, 1427, 1170, 1824, 6417, 1075, 6872, 11]
    sub_seq_size = 150
    decoded_seq = ""
    for i in range(len(hash_keys)):
        sub_seq = sequence[sub_seq_size*i:sub_seq_size*(i+1)]
        decoded_sub_seq = source_decoding.remove_filter(sub_seq, str(hash_keys[i]))
        decoded_seq += decoded_sub_seq
    print(decoded_seq)
   

# =================== main ======================= #
if __name__ == '__main__':
     
    if len(sys.argv) != 2:
        print("usage : sequence_control.py sequence_path")
        sys.exit(1)

    # _, sequence = dfr.read_single_sequence_fasta(sys.argv[1])
    
    #find_hash_keys(sequence)
    #decode_sequence(sequence)
    
    input_sequences = dfr.read_fasta(sys.argv[1])
    
    for fragment_name, fragment in input_sequences.items():
        print("control",fragment_name)
        sequence_check(fragment, window_size=472, verbose=True)

