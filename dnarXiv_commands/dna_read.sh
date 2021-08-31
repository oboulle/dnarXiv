#!/bin/bash

help_function () {
   echo ""
   echo "Usage: dna_read Cname DI Dname"
   echo -e "\tCname : path to the container"
   echo -e "\tDI : index of the document to read"
   echo -e "\tDname : path to save the document"
   echo ""
   exit 1 # Exit script after printing help
}

check_error_function () { #end the program if the previously called script has returned an error
	if [ ! $? = 0 ]
	then
		echo "error in $1"
		exit 1
	fi
}

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

container_path="${1}"
document_index="${2}"
document_path="${3}"

if [ "$#" != 3 ]
then
	echo "error invalid number of arguments"
	help_function
	exit 1
fi

if [ ! -d "$container_path" ] 
then
    echo "error : container $container_path not found" 
    exit 1
fi

if [ ! -d "$container_path/$document_index" ] 
then
    echo "error : document $document_index of $container_path not found" 
    exit 1
fi

#---------------------------------------------#
######### ====== read document ====== #########
#---------------------------------------------#
cdi_file="$container_path/.cdi"
container_index=$(head -n 1 "$cdi_file")

if (( $container_index > 0 ))
then
	echo "the container is not readable"
	exit 1
fi

if [ $HOME == "/Users/oboulle" ]
then
  project_dir="/Users/oboulle/Documents"
else
  project_dir="/home/oboulle/Documents"
fi

# get parameters from the container options
while read var value; do
    export "$var"="$value"
done < "$container_path/.options"

#simulation frag_length spacer 

# get parameters from the .meta file of the document to read
while read var value; do
    export "$var"="$value"
done < "$container_path/$document_index/.meta"

#type n_frag metadata start_primer stop_primer

#----Molecule Selection----#
start_time=$(date +"%s")

nbr_seq=$(($n_frag * 15))
molecule_selection_script="$project_dir/sequencing_simulation/select_sequences.py" 
selected_mol_path="$container_path/$document_index/select_mol.fasta"
#select molecules from container molecules with the good primers
python3 $molecule_selection_script "$container_path/container_molecules.fasta" $selected_mol_path $start_primer $stop_primer $nbr_seq
check_error_function "molecule selection"
seq_sel_time=$(date +"%s")

#----Sequencing & Basecalling----#

#sequencing + basecalling the selected molecules #TODO
convert_fasta_script="$project_dir/synthesis_simulation/dna_file_reader.py" #script to convert fasta to fastq
sequenced_mol_path="$container_path/$document_index/sequenced_mol.fastq"
python3 $convert_fasta_script $selected_mol_path $sequenced_mol_path
check_error_function "sequencing/basecalling"
seq_bc_time=$(date +"%s")

#----Clustering----#

consensus_frag_dir="$container_path/$document_index/frag"
rm -rf consensus_frag_dir
mkdir $consensus_frag_dir

clustering_script="$project_dir/sequencing_simulation/clustering" #script for the clustering
$clustering_script $sequenced_mol_path $consensus_frag_dir "$spacer" $frag_length
check_error_function "clustering"
clust_time=$(date +"%s")

#----Consensus----#

consensus_script="$project_dir/sequencing_simulation/consensus.py"
consensus_path="$container_path/$document_index/consensus.fasta"
python3 $consensus_script "$consensus_frag_dir" $consensus_path "$spacer" $frag_length
check_error_function "consensus"
consensus_time=$(date +"%s")

#----Channel Decoding----#

channel_decoding_script="$project_dir/channel_code/file_decoder.sh"
decoded_sequences_path="$container_path/$document_index/decoded_sequences.fasta"
$channel_decoding_script $consensus_path $frag_length $decoded_sequences_path "$container_path/$document_index/validity_check.txt"
check_error_function "channel decoding"
channel_time=$(date +"%s")

#----Source Decoding----#

source_decoding_script="$project_dir/source_encoding/source_decoding.py"
python3 $source_decoding_script $decoded_sequences_path $document_path $type $frag_length
check_error_function "source decoding"

echo "Document $document_index of $container_path successfully saved to $document_path !"
echo ""
end_time=$(date +"%s")
echo "dna_read : $(($end_time - $start_time)) s" >> workflow_times.txt
echo "   > sequences_select       : $(($seq_sel_time - $start_time)) s" >> workflow_times.txt
echo "   > sequencing_basecalling : $(($seq_bc_time - $seq_sel_time)) s" >> workflow_times.txt
echo "   > clustering             : $(($clust_time - $seq_bc_time)) s" >> workflow_times.txt
echo "   > consensus              : $(($consensus_time - $clust_time)) s" >> workflow_times.txt
echo "   > channel_decoding       : $(($channel_time - $consensus_time)) s" >> workflow_times.txt
echo "   > source_decoding        : $(($end_time - $channel_time)) s" >> workflow_times.txt

exit 0
