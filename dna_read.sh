#!/bin/bash

#DNA READ - extract a document from a container

help_function () {
   echo ""
   echo "Usage: dna_read [-n_read int] Cname DI Dname"
   echo -e "\t-n_read : number of molecules to read [default = 10*n_frag]"
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
case "$1" in
    -n_read ) n_read=$2 ; container_path="${3}"; document_index="${4}"; document_path="${5}";;
    -h | --help ) help_function ; exit 1;;
    -* ) echo "unknown parameter $1" ; exit 1;;
    * ) container_path="${1}"; document_index="${2}"; document_path="${3}";;
esac

if [ "$#" != 3 ] && [ "$#" != 5 ]
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

if [ ! -d "$container_path"/$document_index ] 
then
    echo "error : document $document_index of $container_path not found" 
    exit 1
fi

int_regex='^[0-9]+$'
if test "$n_read" && ! [[ $n_read =~ $int_regex ]]
then
   echo "error: sequencing number ($n_read) is not a number" ; exit 1
fi

#---------------------------------------------#
######### ====== read document ====== #########
#---------------------------------------------#
cdi_file="$container_path"/.cdi
container_index=$(head -n 1 "$cdi_file")

if (( $container_index > 0 ))
then
	echo "the container is not readable"
	exit 1
fi

if [ $HOME == "/Users/oboulle" ]
then
	project_dir="/Users/oboulle/Documents"
elif [ $HOME == "/udd/oboulle" ]
then
	project_dir="/udd/oboulle/Documents"
else
	project_dir="/home/oboulle/Documents"
fi

# get parameters from the container options
while read var value; do
    export "$var"="$value"
done < "$container_path"/.options

#simulation frag_length spacer 

# get parameters from the .meta file of the document to read
stored_document_path="$container_path"/$document_index
while read var value; do
    export "$var"="$value"
done < "$stored_document_path"/.meta

#channel_coding type n_frag start_sequence metadata start_primer stop_primer

#----Molecule Selection----#
start_time=$(date +"%s")

if test -z "$n_read"
then
	#default value of molecule sequencing number
	n_read=$(($n_frag * 50))
fi

molecule_selection_script="$project_dir"/sequencing_simulation/select_molecules.py
selected_mol_path="$stored_document_path"/6_select_mol.fasta
#select molecules from container molecules with the good primers
python3 "$molecule_selection_script" -i "$container_path"/container_molecules.fasta -o "$selected_mol_path" --start $start_primer --stop $stop_primer -n $n_read
check_error_function "molecule selection"
mol_sel_time=$(date +"%s")

#----Sequencing & Basecalling----#

convert_fasta_script="$project_dir"/synthesis_simulation/dna_file_reader.py #script to convert fasta to fastq
simu_seq_bc_script="$project_dir"/sequencing_simulation/sequencing_basecalling_simulator/sequencing_basecalling_simulator.jl
sequenced_reads_path="$stored_document_path"/7_sequenced_reads.fastq
#python3 "$convert_fasta_script" "$selected_mol_path" "$sequenced_reads_path"
"$simu_seq_bc_script" "$selected_mol_path" "$sequenced_reads_path"

check_error_function "sequencing/basecalling"
seq_bc_time=$(date +"%s")

#----Consensus----#

consensus_script="$project_dir"/sequencing_simulation/kmer_consensus/kmer_consensus.py
consensus_path="$stored_document_path"/8_consensus.fasta
expected_length=$(($n_frag * $frag_length))

python3 "$consensus_script" -i "$sequenced_reads_path" -o "$consensus_path" -s "$start_sequence" -e "$expected_length"
check_error_function "consensus"
consensus_time=$(date +"%s")

#----Channel Decoding----#

channel_decoding_script="$project_dir"/channel_code/file_decoder.sh
decoded_sequence_path="$stored_document_path"/9_decoded_sequence.fasta

if $channel_coding
then
	"$channel_decoding_script" "$consensus_path" $frag_length "$decoded_sequence_path" "$container_path"/$document_index/validity_check.txt
else
	cp "$consensus_path" "$decoded_sequence_path"
fi
check_error_function "channel decoding"

channel_time=$(date +"%s")

#----Source Decoding----#

source_decoding_script="$project_dir"/source_encoding/source_decoding.py
reconstructed_source_path="$stored_document_path"/10_reconstructed_source.fasta
python3 "$source_decoding_script" -i "$decoded_sequence_path" -r "$reconstructed_source_path" -o "$document_path" -t $type -m "$metadata"
check_error_function "source decoding"

echo "Document $document_index of $container_path successfully saved to $document_path !"

end_time=$(date +"%s")
times_file="$stored_document_path"/workflow_times.txt
echo "dna_read : $(($end_time - $start_time)) s" >> "$times_file"
echo "   > molecules_select       : $(($mol_sel_time - $start_time)) s" >> "$times_file"
echo "   > sequencing_basecalling : $(($seq_bc_time - $mol_sel_time)) s" >> "$times_file"
echo "   > consensus              : $(($consensus_time - $seq_bc_time)) s" >> "$times_file"
echo "   > channel_decoding       : $(($channel_time - $consensus_time)) s" >> "$times_file"
echo "   > source_decoding        : $(($end_time - $channel_time)) s" >> "$times_file"

exit 0
