
import os
import sys
import inspect
from datetime import datetime

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, os.path.dirname(parentdir)+"/synthesis_simulation")
import dna_file_reader as dfr


def get_size_and_precision(stored_document_path):
    """
    return the size of the document to encode and the result precision percentage with 2 digit after comma
    """
    
    source_path = stored_document_path+"/0_source.fasta"
    result_path = stored_document_path+"/11_reconstructed_source.fasta"
    
    if not os.path.isfile(source_path):
        return "None", "None"
    
    _, source = dfr.read_single_sequence_fasta(source_path)

    if not os.path.isfile(result_path):
        return len(source), "None"
      
    _, result = dfr.read_single_sequence_fasta(result_path)
    
    source_length = len(source)
    nbr_match = 0 #number of matching nucleotides between source and result
    for i in range(source_length):
        if source[i] == result[i]:
            nbr_match += 1
    return len(source), round(100*nbr_match/source_length, 2)

def get_mol_number(stored_document_path):
    """
    number of molecules sequenced in the reading phase
    """
    mol_file_path = stored_document_path+"/6_select_mol.fasta"
    if not os.path.isfile(mol_file_path):
        return "None"
    else:
        return int(sum(1 for line in open(mol_file_path))/2) #number of lines divided by 2
    

def get_reading_time(stored_document_path):
    times_file_path = stored_document_path+"/workflow_times.txt"

    if not os.path.isfile(times_file_path):
        return "None"
    
    times_file = open(times_file_path)
    line = times_file.readline()
    reading_time = "None"
    while line != "":
        if line.startswith("dna_read"):
            reading_time = line.replace("dna_read : ", "").replace(" s\n", "")
            break
        line = times_file.readline()
    times_file.close()
    return reading_time

def write_results_to_file(output_path):
    output = open(output_path, "a") #append at the end of the document
    output.write("date "+current_date+"\n")
    output.write("source_length "+str(source_size)+"\n")
    output.write("n_mol "+str(n_mol)+"\n")
    output.write("precision "+str(precision)+"\n")
    output.write("reading_time "+str(sum_times)+"\n")
    output.write("___\n")
    output.close()

# =================== main ======================= #
if __name__ == '__main__':
    
    if len(sys.argv) != 3:
        print("usage : result_analysis.py stored_document_path output_path")
        sys.exit(1)

    stored_document_path = sys.argv[1]
    
    if not os.path.isdir(stored_document_path):
        print("error :",stored_document_path,"is not a directory")
        sys.exit(1)
        
    print("result analysis...")
 
    source_size, precision = get_size_and_precision(stored_document_path)
    n_mol = get_mol_number(stored_document_path)
    sum_times = get_reading_time(stored_document_path)
    
    current_date = datetime.now().strftime("%d/%m/%Y_%H:%M:%S")
    
    output_path = sys.argv[2]

    write_results_to_file(output_path)
    
    print("\tcompleted !")
