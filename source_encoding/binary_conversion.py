import sys
import hashing

"""
Can be used to convert any type of file
use the byte encoding of the file to convert to a dna sequence
a single substitution error in the sequences can lead to a corrupted file (except for .txt files)
"""

# pair of bit to nucleotide
bit_pair_to_dna = {"00": "A", "01": "G", "10": "T", "11": "C"}
# "_" is for missing nucleotides
dna_to_bit_pair = {"A": "00", "G": "01", "T" :"10", "C": "11", "_" :"00"}

bit_to_dna_x = {"0": "A", "1": "T"}
bit_to_dna_y = {"0": "G", "1": "C"}

dna_to_bit_xy = {"A": "0", "T": "1", "G": "0",  "C": "1", "_": "0"}


def binary_to_dna(binary_string: str) -> str:
    """
    basic conversion of binaries into dna sequence
    conversion is "00"->"A", "01"->"G", "10"->"T", "11"->"C"
    """        
    sequence = ""
    for i in range(0, len(binary_string), 2):
        bit2 = binary_string[i:i+2]
        nucleotide = bit_pair_to_dna[bit2]
        sequence += nucleotide
    return sequence


def dna_to_binary(sequence: str) -> str:
    """
    convert back a dna sequence into binaries
    """
    binary_string = ""
    for nucleotide in sequence:
        bit2 = dna_to_bit_pair[nucleotide]
        binary_string += bit2
    return binary_string


def binary_to_dna_xy(binary_string: str) -> str:
    """
    convert binaries into dna sequence with some properties
    the binaries are divided into parts of 6 bits
    for each 6 bits part, 1st bit : (0:A, 1:T)
                          2nd bit : (0:G, 1:C)
                          3rd 4th bits are normally converted into dna : (00:A, 01:G, 10:T, 11:C)
                          5th 6th bits are normally converted into dna
    this ensure that the total sequence cannot contain homopolymeres > 3,
    the GC% is in [25%, 75%] and the conversion rate is 1.5bit/base
    """
    
    sequence = ""
    n_sextuplet, rest = divmod(len(binary_string), 6)
    for i in range(0, n_sextuplet):
        bits6 = binary_string[i*6:(i+1)*6]
        nucleotide_1 = bit_to_dna_x[bits6[0]]
        nucleotide_2 = bit_to_dna_y[bits6[1]]
        nucleotide_3 = bit_pair_to_dna[bits6[2:4]]
        nucleotide_4 = bit_pair_to_dna[bits6[4:6]]    
        sequence += nucleotide_1 + nucleotide_2 + nucleotide_3 + nucleotide_4
        
    bit_rest = binary_string[n_sextuplet*6:n_sextuplet*6+rest] # rest is 0-2-4 bits
    if len(bit_rest) >= 2: 
        sequence += bit_to_dna_x[bit_rest[0]] + bit_to_dna_y[bit_rest[1]]
        if len(bit_rest) == 4:
            sequence += bit_pair_to_dna[bit_rest[2:4]]
    
    return sequence


def dna_to_binary_xy(sequence: str) -> str:
    """
    convert back a sequence to binaries (opposite of binary_to_dna_xy)
    #TODO possibility to detect and solve some errors
    """
    
    binary_string = ""
    n_quadruplet, rest = divmod(len(sequence), 4)
    for i in range(0, n_quadruplet):
        nucleotides4 = sequence[i*4:(i+1)*4]
        bit_1 = dna_to_bit_xy[nucleotides4[0]] # should only be A or T #TODO CG equals error
        bit_2 = dna_to_bit_xy[nucleotides4[1]] # should only be G or C #TODO AT equals error
        bits_3456 = dna_to_bit_pair[nucleotides4[2]] + dna_to_bit_pair[nucleotides4[3]]
        binary_string += bit_1 + bit_2 + bits_3456
    
    sequence_rest = sequence[n_quadruplet*4:n_quadruplet*4+rest]
    if len(sequence_rest) >= 1:
        binary_string += dna_to_bit_xy[sequence_rest[0]]
        if len(sequence_rest) >= 2:
            binary_string += dna_to_bit_xy[sequence_rest[1]]
            if len(sequence_rest) >= 3:
                binary_string += dna_to_bit_pair[sequence_rest[2]]
    
    return binary_string


def binary_to_dna_z(binary_string: str) -> str:
    """
    convert binaries into dna sequence with some properties
    the binaries are divided into parts of 8 bits
    for each 8 bits part, 1st 2nd bits are normally converted into dna : (00:A, 01:G, 10:T, 11:C)
                          3rd bit is converted depending on previous conversion : if previous in {AT}, then (0:G, 1:C), elif in {GC}, then (0:A, 1:T)
                          4th 5th bits are normally converted into dna
                          6th 7th bits are normally converted into dna
                          8th bit is converted depending on previous conversion : if previous in {AT}, then (0:G, 1:C), elif in {GC}, then (0:A, 1:T)
    this ensure that the total sequence cannot contain homopolymeres > 3,
    the GC% is in [40%, 60%] and the conversion rate is 1.6bit/base
    """
    
    sequence = ""
    n_octets, rest = divmod(len(binary_string), 8)
    for i in range(0, n_octets):
        bits8 = binary_string[i*8:(i+1)*8]
        nucleotide_1 = bit_pair_to_dna[bits8[0:2]]
        if nucleotide_1 in ["A", "T"]:
            nucleotide_2 = bit_to_dna_y[bits8[2]]
        else:
            nucleotide_2 = bit_to_dna_x[bits8[2]]
        nucleotide_3 = bit_pair_to_dna[bits8[3:5]]
        nucleotide_4 = bit_pair_to_dna[bits8[5:7]]
        if nucleotide_4 in ["A", "T"]:
            nucleotide_5 = bit_to_dna_y[bits8[7]]
        else:
            nucleotide_5 = bit_to_dna_x[bits8[7]]
        sequence += nucleotide_1 + nucleotide_2 + nucleotide_3 + nucleotide_4 + nucleotide_5
    
    # rest should be 0 because all documents contains a round number of octet
    # but some "0" can be added to fill the fragments
    bit_rest = binary_string[n_octets*8:n_octets*8+rest] # rest is 0-2-4-6 bits
    
    # last conversions depends on the length of the rest
    if len(bit_rest) == 1:
        sequence += bit_to_dna_x[bit_rest[0]] # nucleotide_1 from a single bit
    elif len(bit_rest) >= 2:
        nucleotide_1 = bit_pair_to_dna[bit_rest[0:2]] # nucleotide_1 from a pair of bits
        sequence += nucleotide_1
        if len(bit_rest) >= 3:
            if nucleotide_1 in ["A", "T"]:
                sequence += bit_to_dna_y[bits8[2]] # nucleotide_2 from a single bit and depending on nucleotide 1
            else:
                sequence += bit_to_dna_x[bits8[2]] # nucleotide_2 from a single bit and depending on nucleotide 1
            
            if len(bit_rest) == 4:
                sequence += bit_to_dna_x[bit_rest[3]] # nucleotide_3 from a single bit
            elif len(bit_rest) >= 5:
                sequence += bit_pair_to_dna[bit_rest[3:5]] # nucleotide_3 from a pair of bits
                
                if len(bit_rest) == 6:
                    sequence += bit_to_dna_x[bit_rest[5]] # nucleotide_4 from a single bit
                elif len(bit_rest) == 7:
                    sequence += bit_pair_to_dna[bit_rest[5:7]] # nucleotide_4 from a pair of bits
    
    return sequence
    

def dna_to_binary_z(sequence: str) -> str:
    """
    convert back a sequence to binaries (opposite of binary_to_dna_z)
    #TODO possibility to detect and solve some errors (not possible if using the ban words removing algorithm in synthesis_modules.ordered_fragmentation)
    """
    binary_string = ""
    n_quintuplet, rest = divmod(len(sequence), 5)
    for i in range(0, n_quintuplet):
        nucleotides5 = sequence[i*5:(i+1)*5]
        bits_12 = dna_to_bit_pair[nucleotides5[0]]
        bit_3 = dna_to_bit_xy[nucleotides5[1]]
        bits_45 = dna_to_bit_pair[nucleotides5[2]]
        bits_67 = dna_to_bit_pair[nucleotides5[3]]
        bit_8 = dna_to_bit_xy[nucleotides5[4]]
        binary_string += bits_12 + bit_3 + bits_45 + bits_67 + bit_8
    
    sequence_rest = sequence[n_quintuplet*5:n_quintuplet*5+rest]
    if len(sequence_rest) >= 1:
        binary_string += dna_to_bit_pair[sequence_rest[0]]
        if len(sequence_rest) >= 2:
            binary_string += dna_to_bit_xy[sequence_rest[1]]
            if len(sequence_rest) >= 3:
                binary_string += dna_to_bit_pair[sequence_rest[2]]
                if len(sequence_rest) >= 4:
                    binary_string += dna_to_bit_pair[sequence_rest[3]]
    
    return binary_string
    
                    
def apply_binary_filter(binary_string: str) -> str:
    """
    apply a filter to the binary string
    use "0" as a default hash key, could be any key, just use the same for the un-hashing
    """
    filter = hashing.hash_string_to_formated_base2("0", len(binary_string))

    filtered_binary_string = ""
    for i in range(len(binary_string)):
        filter_bit = int(filter[i])
        filtered_binary_string += str((int(binary_string[i])+filter_bit) % 2)
         
    return filtered_binary_string
  
    
def encode_file(input_path: str, fragment_length: int) -> str:
    """
    convert any type of file into a dna_sequence
    returns the resulting dna sequence
    """
    # read the lines as bytes
    with open(input_path, "rb") as f:
        byte_list = []
        byte = f.read(1)
        while byte:
            byte_list.append(bin(ord(byte)))
            byte = f.read(1)
    # convert the bytes into binaries
    binary_string = ""
    for x in byte_list:
        binary_char = x[2:]
        #add some 0 at the beginning of binaries shorter than 8
        while len(binary_char) < 8:
            binary_char = "0"+binary_char
        binary_string+=binary_char
    
    # add some "0" at the end of binary string to get a round number of fragments
    conversion_rate = 1.6 # binary_to_dna_z gives 1.6 bit/base
    overhang_length = 0
    primer_length = 20
    
    fragments_bits = (fragment_length - overhang_length) * conversion_rate # number of bits per fragments, remove overhang length from fragment length
    first_fragment_bits = (fragment_length - primer_length - overhang_length) * conversion_rate # number of bits for the first fragment (has a primer and an overhang)
    last_fragment_bits = (fragment_length - primer_length) * conversion_rate # last forward fragment has a primer but no overhang
    
    rest = (len(binary_string) - first_fragment_bits - last_fragment_bits) % (fragments_bits) # rest of bits
    if rest > 0: 
        binary_string += "0" * int(fragments_bits - rest) # fill with 0 to have a round number of fragments
    
    # apply a filter to the binary string -> shuffle the data to avoid long rows of 0 or 1, and avoid rows repetitions 
    filtered_binary_string = apply_binary_filter(binary_string)
    
    # convert binaries into dna sequence
    sequence = binary_to_dna_z(filtered_binary_string)
        
    return sequence


def decode_file(sequence: str, output_path: str) -> None:
    """
    convert a dna_sequence into the original file, assuming no errors in the sequence
    """
    # convert a dna sequence into binaries
    binary_string = dna_to_binary_z(sequence)
    
    # apply the same filter used in the encoding to the binary string to remove it  
    filtered_binary_string = apply_binary_filter(binary_string)

    # case binaries length is not multiple of 8 -> remove the excess bits
    rest = len(filtered_binary_string) % 8
    if rest != 0:
        filtered_binary_string = filtered_binary_string[:-rest]

    #remove zeros at the end (the last fragment is usually filled with A = 8 zeros
    while filtered_binary_string.endswith(8*"0"):
        filtered_binary_string = filtered_binary_string[:-8]       

    if not filtered_binary_string:
        print("warning file conversion, decoding an empty file")
        return
    # convert binaries into bytes
    n = int(filtered_binary_string, 2)
    bytes = n.to_bytes((n.bit_length() + 7) // 8, 'big')
    decoded_bytes = bytes.decode("utf-8", "ignore")
    # write the bytes into the file
    f = open(output_path, "w")
    f.write(decoded_bytes)
    f.close()


# =================== main ======================= #
if __name__ == '__main__':
    #doc_path = sys.argv[1]
    binary_string = sys.argv[1]
    #seq = encode_file(doc_path)
    seq = binary_to_dna_z(binary_string)
    print(seq)
    print(dna_to_binary_z(seq))
    #decode_file(seq, "test")

