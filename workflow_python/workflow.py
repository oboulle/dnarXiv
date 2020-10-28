import sys
import os
import shutil
from shutil import copyfile

"""
    get the paths of the required scripts and directory from a file
"""
def get_paths(paths_file):
    ff = open(paths_file)
    paths_dict = {}       
    line = ff.readline().replace("\n","")
    while line != "":
        ll = line.split("=")
        paths_dict[ll[0]] = ll[1]
        line = ff.readline().replace("\n","")
    ff.close()
    return paths_dict

"""
    get input parameters from the user to create random sequences
"""
def generate_sequences(base_seq_path):
    nbr_sequence = input(" > nombre de séquences : ")
    size_sequence = input(" > taille des séquences : ")
    h_max = input(" > taille maximale des homopolymères : ")
    return os.system("python3 "+paths_dict["sequence_generator"]+" "+base_seq_path+" "+nbr_sequence+" "+size_sequence+" "+h_max)
    
"""
    get the name of the base sequence file from the user, repeat until valid name is given
"""
def get_seq_file(base_seq_path):
    seq_path = input("chemin du fichier .fasta des séquences : ")
    if not os.path.isfile(seq_path):
        print ("Le fichier est introuvable : "+seq_path)
        get_seq_file()
    else:
        copyfile(seq_path, base_seq_path)
     
"""
    add a primer at the beginning and the end of the sequences
"""
def add_primer(base_seq_path, primer_seq_path):
    sequence_primer_dir = paths_dict["sequence_primer_dir"]
    return os.system("python3 "+sequence_primer_dir+"sequence_primer.py "+base_seq_path+" "+ \
                     primer_seq_path+" "+sequence_primer_dir+"test_primer.fasta")

"""
    get input parameters from the user to simulate the synthesis of the sequences
"""
def synthetise(primer_seq_path, synthesis_path):
    nbr_synth = input(" > nombre de synthèse pour chaque séquence : ")
    i_error = input(" > taux d'erreur d'insertion : ")
    d_error = input(" > taux d'erreur de deletion : ")
    s_error = input(" > taux d'erreur de substitution : ")
    return os.system("python3 "+paths_dict["synthesis_simulation"]+" -i "+primer_seq_path+" -o "+synthesis_path+ \
                     " -n "+nbr_synth+" --i_error "+i_error+" --d_error "+d_error+" --s_error "+s_error)

"""
    sequence the data from the synthesis
"""
def sequencing(synthesis_path, sequencing_path):
    nbr_seq = input(" > nombre de séquençage : ")
    return os.system(paths_dict["deep_simulator"]+" -i "+synthesis_path+" -o "+sequencing_path+" -G 1 -H "+paths_dict["DeepSimulator"]+" -n "+nbr_seq)
  
#________________début du processus________________#

print("___Début du processus___")

paths_dict = get_paths("/udd/oboulle/Documents/workflow_global/workflow_python/script_paths.txt")
process_name = input(" nom du processus : ")
try:
    os.mkdir(process_name)
except OSError:
    shutil.rmtree(process_name)
    os.mkdir(process_name)
    
#________________initialisation des séquences de base________________#

input("___Initialisation des séquences de base___")
base_seq_path = process_name+"/1_base_seq_file.fasta"
random_seq_input = input(" générer des séquences aléatoires ? (y/n) : ")

if random_seq_input == "y" or random_seq_input == "yes":
    result = generate_sequences(base_seq_path)
    if(result != 0):
        sys.exit(1)
    else:
        print("\n séquences générées !\n")
else:
    get_seq_file(base_seq_path)
    
#________________ajout des primers________________#

input("___Ajout de primers en début et fin de séquences___")
primer_seq_path = process_name+"/2_primer_seq_file.fasta"
result = add_primer(base_seq_path, primer_seq_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n primers ajoutés !\n")

#________________synthèse________________#

input("___Synthèse des séquences___")
synthesis_path = process_name+"/3_synthesis_file.fasta"
result = synthetise(primer_seq_path, synthesis_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n synthèse effectuée !\n")

#________________séquençage________________#

input("___Séquençage___")
sequencing_path = process_name+"/4_sequencing"
result = sequencing(synthesis_path, sequencing_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n séquençage effectuée !\n")


print("___Fin du processus !___")







































