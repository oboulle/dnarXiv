import os
import sys
import inspect
from theano.gpuarray.dnn import get_precision

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

def sum_times(times_file_path):
    
    return 0

# =================== main ======================= #
if __name__ == '__main__':
     
    if len(sys.argv) != 6:
        print("usage : result_analysis.py source_path result_path n_seq times_file_path output_path")
        sys.exit(1)

    print("result analysis...")

    source_path = sys.argv[1]
    result_path = sys.argv[2]
    n_seq = int(sys.argv[3])
    times_file_path = sys.argv[4]
    output_path = sys.argv[5]

    _, source = dfr.read_single_sequence_fasta(source_path)
    _, result = dfr.read_single_sequence_fasta(result_path)
    
    precision = get_precision(source, result)
    sum_times = sum_times(times_file_path)
    
    output = open(output_path, "a")
    output.write("source_length "+str(len(source))+"\n")
    output.write("n_seq "+str(n_seq)+"\n")
    output.write("precision "+str(precision)+"\n")
    output.write("time "+str(sum_times)+"\n")
    output.close()
    
    print("\tcompleted !")
