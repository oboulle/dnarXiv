#!/usr/bin/python3

import random
import argparse
import numpy
import dna_file_reader as dfr


"""
simulations of fragment assembly into molecules
"""

def unordered_fragment_assembly(input_path: str, output_path: str, spacer: str, start_primer: str, stop_primer: str, nbr_mol: int) -> None:
    """
    create molecules consisting of the start primer, multiple random fragments separated with the spacer and the stop primer
    the resulting molecules are saved to a .fasta file to simulate the storing
    """
    input_sequences_list = list(dfr.read_fasta(input_path).items())
    file_output = open(output_path, "w+")
    
    mean_n_frag = 10 #TODO set as optional param
    seq_tot = 0
    for i in range(nbr_mol):
        sequence = start_primer
        
        n_frag = max(1,int(numpy.random.normal(mean_n_frag, mean_n_frag/5))) #normal distribution of the number of fragment in generated molecules 
        for j in range(n_frag):
            name, value = random.choice(input_sequences_list)
            sequence += spacer+value
        sequence += spacer+stop_primer
        file_output.write(">molecule_" + str(seq_tot) + "\n")
        file_output.write(sequence + "\n")
        seq_tot += 1
    file_output.close()
    
    
def ordered_fragment_assembly(input_path: str, output_path: str, nbr_fragments: int, nbr_mol: int) -> None:
    """
    simulate an ordered assembly of fragments surrounded by the primers
    input : file of synthesized fragments, names are ">frag_X_Y" with X the fragment number
    the resulting molecules are saved to a .fasta file to simulate the storing
    """
    
    input_sequences_dict = dfr.read_fasta(input_path)
    if len(input_sequences_dict) == 0:
        print("error : not fragments found at :",input_path)
        exit(1)


    frag_sequences_dict = {} # dict of a list of synthesized sequences for each frag number
    for i in range(nbr_fragments):
        # get the list of every synthesized sequence corresponding to the same fragment
        next_frag_list_Fw = [value for key, value in input_sequences_dict.items() if key.startswith("frag_"+str(i)+"_Fw")]
        next_frag_list_Rv = [value for key, value in input_sequences_dict.items() if key.startswith("frag_"+str(i)+"_Rv")]

        frag_sequences_dict[str(i)+"_Fw"] = next_frag_list_Fw
        frag_sequences_dict[str(i)+"_Rv"] = next_frag_list_Rv

    designed_sequences = []
    mol_count = 0
    length_dict = {} # count the frag number of added mols
    assembly_rate = 0.8 # 80% chance to be followed by the next fragment

    while mol_count < nbr_mol:
        new_mol = "" # molecule being assembled
        frag_count = 0
        
        # 1/2 chance to assemble a reverse molecule
        if random.random() < .5:
            direction = "Fw"
        else: 
            direction = "Rv"
        
        for i in range(nbr_fragments):
            if direction == "Fw":
                fragment = random.choice(frag_sequences_dict[str(i)+"_Fw"])
            else:
                fragment = random.choice(frag_sequences_dict[str(nbr_fragments-1-i)+"_Rv"])
            
            
            # use the assembly rate to determine if the new frag will be sticked to the building molecule
            if random.random() < assembly_rate: 
                new_mol += fragment
                frag_count += 1
            else:
                # the fragment is not sticked to the building molecule
                # add the previous builded molecule
                if nbr_fragments < 4 or (frag_count > 1 and frag_count < nbr_fragments): # more realistic : ignore 1 fragment mols and max fragments mols, doesn't make sense when very low frag_number
                    designed_sequences.append(new_mol)
                    length_dict[frag_count] = length_dict.get(frag_count, 0) + 1
                    mol_count += 1
                # a new molecule is created, starting with this frag
                new_mol = fragment
                frag_count = 1
        
        if nbr_fragments < 4 or (frag_count > 1 and frag_count < nbr_fragments):
            designed_sequences.append(new_mol)
            length_dict[frag_count] = length_dict.get(frag_count, 0) + 1
            mol_count += 1         

    #print(sorted(length_dict.items(), key=lambda x: x[0]))
    
    file_output = open(output_path, "w+")
    for i, sequence in enumerate(designed_sequences):
        file_output.write(">molecule_" + str(i) + "\n")
        file_output.write(sequence + "\n")
    file_output.close()


# =================== main =======================#
if __name__ == '__main__':
    # --------- argument part -------------#
    parser = argparse.ArgumentParser(description='generate molecules from fragments')
    parser.add_argument('-i', action='store', dest='input_path', required=True,
                        help='the input fasta file')
    parser.add_argument('-o', action='store', dest='output_path', required=True,
                        help='the output fasta file')
    parser.add_argument('-t', action='store', dest='assembly_type', required=True,
                    help='assembly type of the encoded document (indexed/ordered)')
    parser.add_argument('-f', action='store', dest='nbr_frag', required=True,
                        type=int, help='number of fragments')
    parser.add_argument('-n', action='store', dest='nbr_mol', required=True,
                        type=int, help='number of molecules')

    # ---------- input list ---------------#
    arg = parser.parse_args()
    
    print("fragment assembly...")
    if arg.assembly_type == "indexed":
        unordered_fragment_assembly(arg.input_path, arg.output_path, arg.spacer, arg.start_primer, arg.stop_primer, arg.nbr_mol)
    else:
        ordered_fragment_assembly(arg.input_path, arg.output_path, arg.nbr_frag, arg.nbr_mol)
    print("\tcompleted !")

