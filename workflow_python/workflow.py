import sys
import os
import shutil
from shutil import copyfile


def getPaths(paths_file):
    
    ff = open(paths_file)
    paths_dict = {}
        
    line = ff.readline().replace("\n","")
    while line != "":
        ll = line.split("=")
        paths_dict[ll[0]] = ll[1]
        line = ff.readline().replace("\n","")
    ff.close()
    return paths_dict

def generateSequences(base_seq_path):
    nbr_sequence = input(" > nombre de séquences : ")
    size_sequence = input(" > taille des séquences : ")
    h_max = input(" > taille maximale des homopolymères : ")
    return os.system("python3 "+paths_dict["sequence_generator"]+" "+base_seq_path+" "+nbr_sequence+" "+size_sequence+" "+h_max)
    
def getSeqFile(base_seq_path):
    seq_path = input("chemin du fichier .fasta des séquences : ")
    if not os.path.isfile(seq_path):
        print ("Le fichier est introuvable : "+seq_path)
        getSeqFile()
    else:
        copyfile(seq_path, base_seq_path)
        
def addTag(base_seq_path, tagged_seq_path):
    sequence_tag_dir = paths_dict["sequence_tag_dir"]
    return os.system("python3 "+sequence_tag_dir+"sequence_tag.py "+base_seq_path+" "+ \
                     tagged_seq_path+" "+sequence_tag_dir+"begin_tag.txt "+sequence_tag_dir+"end_tag.txt")

def synthetise(tagged_seq_path, synthesis_path):
    nbr_synth = input(" > nombre de synthèse pour chaque séquence : ")
    i_error = input(" > taux d'erreur d'insertion : ")
    d_error = input(" > taux d'erreur de deletion : ")
    s_error = input(" > taux d'erreur de substitution : ")
    return os.system("python3 "+paths_dict["synthesis_simulation"]+" -i "+tagged_seq_path+" -o "+synthesis_path+ \
                     " -n "+nbr_synth+" --i_error "+i_error+" --d_error "+d_error+" --s_error "+s_error)

def sequencing(synthesis_path, sequencing_path):
    nbr_seq = input(" > nombre de séquençage : ")
    return os.system(paths_dict["deep_simulator"]+" -i "+synthesis_path+" -o "+sequencing_path+" -G 1 -B 3 -H "+paths_dict["DeepSimulator"]+" -n "+nbr_seq)
    
paths_dict = getPaths("/udd/oboulle/Documents/workflow_global/workflow_python/script_paths.txt")

print("___Début du processus___")

process_name = input(" nom du processus : ")
print("\n")
try:
    os.mkdir(process_name)
except OSError:
    shutil.rmtree(process_name)
    os.mkdir(process_name)
    
#________________initialisation des séquences de base________________#

input("___Initialisation des séquences de base___")
base_seq_path = process_name+"/base_seq_file.fasta"
random_seq_input = input(" générer des séquences aléatoires ? (y/n) : ")

if random_seq_input == "y" or random_seq_input == "yes":
    result = generateSequences(base_seq_path)
    if(result != 0):
        sys.exit(1)
    else:
        print("\n séquences générées !\n")
else:
    getSeqFile(base_seq_path)
    
#________________ajout des tags________________#

input("___Ajout de tags en début et fin de séquences___")
tagged_seq_path = process_name+"/tagged_seq_file.fasta"
result = addTag(base_seq_path, tagged_seq_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n tags ajoutés !\n")

#________________synthèse________________#

input("___Synthèse des séquences___")
synthesis_path = process_name+"/synthesis_file.fasta"
result = synthetise(tagged_seq_path, synthesis_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n synthèse effectuée !\n")

#________________séquençage________________#

input("___Séquençage___")
sequencing_path = process_name+"/sequencing"
result = sequencing(synthesis_path, sequencing_path)
if(result != 0):
    sys.exit(1)
else:
    print("\n séquençage effectuée !\n")


print("___Fin du processus !___")







































