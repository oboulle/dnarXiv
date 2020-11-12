import sys
import os
import shutil
from shutil import copyfile

global size_sequence

def get_paths(working_dir, paths_file):
    """
        get the paths of the required projects from a file
    """
    ff = open(paths_file)
    paths_dict = {}       
    line = ff.readline().replace("\n","")
    while line != "":
        ll = line.split("=")
        paths_dict[ll[0]] = working_dir+ll[1]
        line = ff.readline().replace("\n","")
    ff.close()
    return paths_dict

def generate_sequences(base_seq_path):
    """
        get input parameters from the user to create random sequences
    """
    nbr_sequence = input(" > nombre de séquences : ")
    global size_sequence
    size_sequence = input(" > taille des séquences : ")
    h_max = input(" > taille maximale des homopolymères : ")
    return os.system("python3 "+paths_dict["sequence_generator"]+"sequence_generator.py '"+base_seq_path+"' "+nbr_sequence+" "+size_sequence+" "+h_max)

def get_seq_file(base_seq_path):
    """
        get the name of the base sequence file from the user, repeat until valid name is given
    """
    seq_path = input("chemin du fichier .fasta des séquences : ")
    if not os.path.isfile(seq_path):
        print ("Le fichier est introuvable : "+seq_path)
        get_seq_file()
    else:
        copyfile(seq_path, base_seq_path)

def add_primer(base_seq_path, primer_seq_path, process_path):
    """
        add a primer at the beginning and the end of the sequences
    """
    return os.system("python3 "+paths_dict["sequence_primer"]+"sequence_primer.py '"+base_seq_path+"' '"+ \
                     primer_seq_path+"' "+process_path)

def synthetise(primer_seq_path, synthesis_path):
    """
        get input parameters from the user to simulate the synthesis of the sequences
    """
    nbr_synth = input(" > nombre de synthèse pour chaque séquence : ")
    i_error = input(" > taux d'erreur d'insertion : ")
    d_error = input(" > taux d'erreur de deletion : ")
    s_error = input(" > taux d'erreur de substitution : ")
    return os.system("python3 "+paths_dict["synthesis_simulation"]+"synthesis_simulator.py -i '"+primer_seq_path+"' -o '"+synthesis_path+ \
                     "' -n "+nbr_synth+" --i_error "+i_error+" --d_error "+d_error+" --s_error "+s_error)

def sequencing(synthesis_path, sequencing_path):
    """
        sequence the data from the synthesis
    """
    nbr_read = input(" > nombre de lecture : ")
    return os.system(paths_dict["deep_simulator"]+"deep_simulator.sh -i '"+synthesis_path+"' -o '"+sequencing_path+"' -H "+paths_dict["deep_simulator"]+" -C '"+conda_env+"' -n "+nbr_read)

def basecalling(sequencing_path, basecalling_path):
    """
        apply the guppy basecalling
    """
    os.mkdir(basecalling_path)
    return os.system(paths_dict["guppy_basecaller"]+" -r --input_path '"+sequencing_path \
                     +"' --save_path '"+basecalling_path+"' -c dna_r9.4.1_450bps_hac.cfg "\
                     +"--cpu_threads_per_caller 8 --num_callers 1")

def demultiplexing(basecalling_path, demultiplexing_path, process_path):
    """
        sort the sequences by primers
    """
    for file in os.listdir(basecalling_path):
        if file.endswith(".fastq"):
            fastq_sequences_path = os.path.join(basecalling_path, file)
            continue
    return os.system("python3 "+paths_dict["demultiplexing"]+"demultiplexing.py -i '"+fastq_sequences_path+"' -o '"+demultiplexing_path+"' -p '"+process_path+"'/primers")

def consensus(demultiplexing_path, consensus_path, process_path):
    """
        create a consensus for the sequences
    """
    os.mkdir(consensus_path)
    for file in os.listdir(demultiplexing_path):
        if file != "unlinked.fastq":
            fastq_sequence_path = os.path.join(demultiplexing_path, file)
            sequence_name = file.replace(".fastq","")
            result = os.system("python3 " + paths_dict["consensus"] + "ccsa.py -read '"+fastq_sequence_path +"' -primer '"+process_path+"'/primers/"+sequence_name+".fasta "\
            +" -length "+size_sequence+" -out '"+consensus_path+"'/"+sequence_name+".fasta -graph '"+consensus_path+"'/"+sequence_name+".gexf")
            if result != 0:
                return result
    return 0

#________________Début du processus________________#
running_path = os.getcwd() #place where this script is run
if "genouest" in running_path:
    working_dir = "/home/genouest/genscale/oboulle/documents/"
    conda_env = "/home/genouest/genscale/oboulle/anaconda2"
else:
    working_dir="/home/oboulle/Documents/"
    conda_env="/home/oboulle/anaconda2"
os.chdir(working_dir+"workflow_global/workflow_1")

print("___Début du processus___")

paths_dict = get_paths(working_dir, "project_paths.txt")
process_path = running_path+"/"+input(" nom du processus : ")

try:
    os.mkdir(process_path)
except OSError:
    shutil.rmtree(process_path)
    os.mkdir(process_path)
    
#________________Initialisation des séquences de base________________#

input("___Initialisation des séquences de base___")
base_seq_path = process_path+"/1_base_seq_file.fasta"
random_seq_input = input(" générer des séquences aléatoires ? (y/n) : ")

if random_seq_input == "y" or random_seq_input == "yes":
    result = generate_sequences(base_seq_path)
    if(result != 0):
        sys.exit(1)
    else:
        print("\n séquences générées !\n")
else:
    get_seq_file(base_seq_path)
    
#________________Ajout des primers________________#

input("___Ajout de primers en début et fin de séquences___")
primer_seq_path = process_path+"/2_primer_seq_file.fasta"
result = add_primer(base_seq_path, primer_seq_path, process_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n primers ajoutés !\n")

#________________Synthèse________________#

input("___Synthèse des séquences___")
synthesis_path = process_path+"/3_synthesis_file.fasta"
result = synthetise(primer_seq_path, synthesis_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n synthèse effectuée !\n")

#________________Séquençage________________#

input("___Séquençage___")
sequencing_path = process_path+"/4_sequencing"
result = sequencing(synthesis_path, sequencing_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n séquençage effectué !\n")
    
#________________Base Calling________________#

input("___Base Calling___")
basecalling_path = process_path+"/5_basecalling"
result = basecalling(sequencing_path, basecalling_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n base calling effectué !\n")

#________________Demultiplexing________________#

input("___Demultiplexing___")
demultiplexing_path = process_path+"/6_demultiplexing"
result = demultiplexing(basecalling_path, demultiplexing_path, process_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n demultiplexing effectué !\n")

#________________Consensus________________#
input("___Consensus___")
consensus_path = process_path+"/7_consensus"
result = consensus(demultiplexing_path, consensus_path, process_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n consensus effectué !\n")

#________________Fin________________#
print("___Fin du processus !___")






