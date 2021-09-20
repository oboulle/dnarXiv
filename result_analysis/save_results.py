
import os
import sys
import inspect
from datetime import datetime

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, os.path.dirname(parentdir)+"/synthesis_simulation")
import dna_file_reader as dfr


def get_precision(source, result):
    """
    return the precision percentage with 2 digit after comma
    """
    source_length = len(source)
    nbr_match = 0 #number of matching nucleotides between source and result
    for i in range(source_length):
        if source[i] == result[i]:
            nbr_match += 1
    return round(100*nbr_match/source_length, 2)

def get_reading_time(times_file_path):
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

    source_path = stored_document_path+"/0_source.fasta"
    result_path = stored_document_path+"/11_reconstructed_source.fasta"
    n_mol = int(sum(1 for line in open(stored_document_path+"/6_select_mol.fasta"))/2) #number of lines divided by 2
    times_file_path = stored_document_path+"/workflow_times.txt"
    output_path = sys.argv[2]

    _, source = dfr.read_single_sequence_fasta(source_path)
    _, result = dfr.read_single_sequence_fasta(result_path)
    
    precision = get_precision(source, result)
    sum_times = get_reading_time(times_file_path)
    
    now = datetime.now()
    current_date = now.strftime("%d/%m/%Y_%H:%M:%S")
    
    output = open(output_path, "a") #append a the end of the document
    output.write("date "+current_date+"\n")
    output.write("source_length "+str(len(source))+"\n")
    output.write("n_mol "+str(n_mol)+"\n")
    output.write("precision "+str(precision)+"\n")
    output.write("reading_time "+str(sum_times)+"\n")
    output.write("___\n")
    output.close()
    
    print("\tcompleted !")
