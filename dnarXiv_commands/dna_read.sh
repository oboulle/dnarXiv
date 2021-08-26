#!/bin/bash

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

help_function() {
   echo ""
   echo "Usage: dna_read Cname DI Dname"
   echo -e "\tCname : path to the container"
   echo -e "\tDI : index of the document to read"
   echo -e "\tDname : path to save the document"
   echo ""
   exit 1 # Exit script after printing help
}

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

#----Sequence Selection----#

nbr_seq=50 #TODO
sequence_selection_script="$project_dir/sequencing_simulation/select_sequences.py" 

#sequence molecules from container molecules with the good primers
python3 $sequence_selection_script "$container_path/container_molecules.fasta" "$container_path/$document_index/select_mol.fasta" $start_primer $stop_primer $nbr_seq

#----Sequencing & Basecalling----#

#sequencing + basecalling the selected molecules #TODO
convert_fasta_script="$project_dir/synthesis_simulation/dna_file_reader.py" #script to convert fasta to fastq
python3 $convert_fasta_script "$container_path/$document_index/select_mol.fasta" "$container_path/$document_index/sequenced_mol.fastq"

#----Clustering----#

consensus_frag_dir="$container_path/$document_index/frag"
rm -rf consensus_frag_dir
mkdir $consensus_frag_dir

clustering_script="$project_dir/sequencing_simulation/clustering" #script for the clustering
$clustering_script "$container_path/$document_index/sequenced_mol.fastq" $consensus_frag_dir "$spacer" $frag_length
if [ ! $? = 0 ]
then
	echo "error in clustering"
	exit 1
fi

#----Consensus----#

consensus_script="$project_dir/sequencing_simulation/consensus.py"
consensus_path="$container_path/$document_index/consensus.fasta"
python3 $consensus_script "$consensus_frag_dir" $consensus_path "$spacer" $frag_length

#----Channel Decoding----#

channel_decoding_script="$project_dir/channel_code/file_decoder.sh"
decoded_sequences_path="$container_path/$document_index/decoded_sequences.fasta"
$channel_decoding_script $consensus_path $frag_length $decoded_sequences_path "$container_path/$document_index/validity_check.txt"

#----Source Decoding----#

source_decoding_script="$project_dir/source_encoding/source_decoding.py"
python3 $source_decoding_script $decoded_sequences_path $document_path $type $frag_length

echo "Document $document_index of $container_path successfully saved to $document_path !"
exit 0
