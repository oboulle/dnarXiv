import sys
import os
import shutil
from shutil import copyfile

sys.path.insert(0, '../../synthesis_simulation/utils')
import dna_file_reader as dfr


def get_paths(working_dir, paths_file):
    """
    get the paths of the required projects from a file
    :param working_dir: directory containing the workflow projects
    :param paths_file: path to the file project_paths.txt containing the path of the used projects
    :return: a dictionary containing the required paths for the workflow
    """
    ff = open(paths_file)
    paths_dict = {}
    line = ff.readline().replace("\n", "")
    while line != "":
        ll = line.split("=")
        paths_dict[ll[0]] = working_dir + ll[1]
        line = ff.readline().replace("\n", "")
    ff.close()
    return paths_dict


def generate_sequences(base_seq_path):
    """
    get input parameters from the user to create random sequences
    :param base_seq_path: path to save the created sequences
    :return: the size of the sequences chosen by the user and the result of the python script
    """
    nbr_sequence = input(" > nombre de séquences : ")
    sequences_size = input(" > taille des séquences : ")
    h_max = input(" > taille maximale des homopolymères : ")
    return sequences_size, os.system("python3 " + paths_dict[
        "sequence_generator"] + "sequence_generator.py '" + base_seq_path + "' " + nbr_sequence + " " + sequences_size + " " + h_max)


def get_seq_file(base_seq_path):
    """
    get the name of the initial sequences file from the user, repeat until valid name is given
    copy the file in the base_seq_path
    :param base_seq_path: path to save the .fasta file
    """
    seq_path = input("chemin du fichier .fasta des séquences : ")
    if not os.path.isfile(seq_path):
        print("Le fichier est introuvable : " + seq_path)
        get_seq_file(base_seq_path)
    else:
        copyfile(seq_path, base_seq_path)
        return dfr.get_sequence_size(seq_path)


def add_primer(base_seq_path, primer_seq_path, process_path):
    """
    generate and add a primer at the beginning and the end of the sequences
    :param base_seq_path: path of the .fasta file of the sequences
    :param primer_seq_path: path to save the primered sequences
    :param process_path: directory to create the primers dir
    :return: the result of the python script
    """
    return os.system("python3 " + paths_dict["sequence_primer"] + "sequence_primer.py '" + base_seq_path + "' '" +
                     primer_seq_path + "' " + process_path)


def synthesise(primer_seq_path, synthesis_path):
    """
    get input parameters from the user to simulate the synthesis of the sequences
    :param primer_seq_path: path of the sequences to synthesis (primers added)
    :param synthesis_path: path to save the synthesised sequences
    :return: the result of the python script
    """
    nbr_synth = input(" > nombre de synthèse pour chaque séquence : ")
    i_error = input(" > taux d'erreur d'insertion : ")
    d_error = input(" > taux d'erreur de deletion : ")
    s_error = input(" > taux d'erreur de substitution : ")
    return os.system("python3 " + paths_dict[
        "synthesis_simulation"] + "synthesis_simulator.py -i '" + primer_seq_path + "' -o '" + synthesis_path +
                     "' -n " + nbr_synth + " --i_error " + i_error + " --d_error " + d_error + " --s_error " + s_error)


def sequencing(synthesis_path, sequencing_path):
    """
    sequence the data from the synthesis
    :param synthesis_path: path of the .fasta file to sequence
    :param sequencing_path: directory to save the result of the sequencing
    :return: the result of the shell script
    """
    nbr_read = input(" > nombre de lecture : ")
    return os.system(paths_dict[
                         "deep_simulator"] + "deep_simulator.sh -i '" + synthesis_path + "' -o '" + sequencing_path + "' -H " +
                     paths_dict["deep_simulator"] + " -C '" + conda_env + "' -n " + nbr_read)


def basecalling(sequencing_path, basecalling_path):
    """
    apply the guppy basecalling
    :param sequencing_path: path of the directory of the sequenced files
    :param basecalling_path: directory to save the result of the basecalling
    :return: the result of the basecalling script
    """
    os.mkdir(basecalling_path)
    return os.system(paths_dict["guppy_basecaller"] + " -r --input_path '" + sequencing_path
                     + "' --save_path '" + basecalling_path + "' -c dna_r9.4.1_450bps_hac.cfg "
                     + "--cpu_threads_per_caller 8 --num_callers 1")


def demultiplexing(basecalling_path, demultiplexing_path, process_path):
    """
    sort the sequences by primers
    :param basecalling_path: directory of the sequences to use
    :param demultiplexing_path: directory to save the result of the demultiplexing
    :param process_path: directory containing the primers directory
    :return: the result of the python script
    """
    fastq_sequences_path = None
    for file in os.listdir(basecalling_path):
        if file.endswith(".fastq"):
            fastq_sequences_path = os.path.join(basecalling_path, file)
            continue
    if fastq_sequences_path is None:  # no .fastq file found
        return 1
    return os.system("python3 " + paths_dict[
        "demultiplexing"] + "demultiplexing.py -i '" + fastq_sequences_path + "' -o '" + demultiplexing_path + "' -p '" + process_path + "'/primers")


def consensus(demultiplexing_path, consensus_path, process_path, sequences_size):
    """
    create a consensus for the sequences
    :param demultiplexing_path: directory of the sequences to use
    :param consensus_path: directory to save the result of the consensus
    :param process_path: directory containing the primers directory
    :param sequences_size: size of the initial sequences (before primer addition)
    :return: the result of the python script
    """
    os.mkdir(consensus_path)
    for file in os.listdir(demultiplexing_path):
        if file != "unlinked.fastq":
            fastq_sequence_path = os.path.join(demultiplexing_path, file)
            sequence_name = file.replace(".fastq", "")
            result = os.system("python3 " + paths_dict[
                "consensus"] + "ccsa.py -read '" + fastq_sequence_path + "' -primer '" + process_path + "'/primers/" + sequence_name + ".fasta "
                               + " -length " + sequences_size + " -out '" + consensus_path + "'/" + sequence_name + ".fasta -graph '" + consensus_path + "'/" + sequence_name + ".gexf")
            if result != 0:
                return result
    return 0


# ________________Start of the process________________ #
running_path = os.getcwd()  # place where this script is run
if "genouest" in running_path:
    working_dir = "/home/genouest/genscale/oboulle/documents/"
    conda_env = "/home/genouest/genscale/oboulle/anaconda2"
else:
    working_dir = "/home/oboulle/Documents/"
    conda_env = "/home/oboulle/anaconda2"
os.chdir(working_dir + "workflow_global/workflow_1")

print("___Début du processus___")

paths_dict = get_paths(working_dir, "project_paths.txt")
process_path = running_path + "/" + input(" nom du processus : ")

try:
    os.mkdir(process_path)
except OSError:
    shutil.rmtree(process_path)
    os.mkdir(process_path)

# ________________Initialisation of base sequences________________ #

input("___Initialisation des séquences de base___")
base_seq_path = process_path + "/1_base_seq_file.fasta"
random_seq_input = input(" générer des séquences aléatoires ? (y/n) : ")

if random_seq_input == "y" or random_seq_input == "yes":
    sequences_size, result = generate_sequences(base_seq_path)
    if result != 0:
        sys.exit(1)
    else:
        print("\n séquences générées !\n")
else:
    sequences_size = get_seq_file(base_seq_path)

# ________________Primers addition________________ #

input("___Ajout de primers en début et fin de séquences___")
primer_seq_path = process_path + "/2_primer_seq_file.fasta"
result = add_primer(base_seq_path, primer_seq_path, process_path)
if result != 0:
    sys.exit(1)
else:
    print("\n primers ajoutés !\n")

# ________________Synthesis________________ #

input("___Synthèse des séquences___")
synthesis_path = process_path + "/3_synthesis_file.fasta"
result = synthesise(primer_seq_path, synthesis_path)
if result != 0:
    sys.exit(1)
else:
    print("\n synthèse effectuée !\n")

# ________________Sequencing________________ #

input("___Séquençage___")
sequencing_path = process_path + "/4_sequencing"
result = sequencing(synthesis_path, sequencing_path)
if result != 0:
    sys.exit(1)
else:
    print("\n séquençage effectué !\n")

# ________________Base Calling________________ #

input("___Base Calling___")
basecalling_path = process_path + "/5_basecalling"
result = basecalling(sequencing_path, basecalling_path)
if result != 0:
    sys.exit(1)
else:
    print("\n base calling effectué !\n")

# ________________Demultiplexing________________ #

input("___Demultiplexing___")
demultiplexing_path = process_path + "/6_demultiplexing"
result = demultiplexing(basecalling_path, demultiplexing_path, process_path)
if result != 0:
    sys.exit(1)
else:
    print("\n demultiplexing effectué !\n")

# ________________Consensus________________ #

input("___Consensus___")
consensus_path = process_path + "/7_consensus"
result = consensus(demultiplexing_path, consensus_path, process_path, sequences_size)
if result != 0:
    sys.exit(1)
else:
    print("\n consensus effectué !\n")

# ________________End________________ #
print("___Fin du processus !___")
