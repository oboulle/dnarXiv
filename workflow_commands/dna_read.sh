#!/bin/bash
set -u #exit and display error message if a variable is empty 

#DNA READ - extract a document from a container

project_dir="$(dirname $(dirname ${BASH_SOURCE}))" #parent of the directory containing this script
source "$project_dir"/workflow_commands/metadata_manager.sh #load the xml manager script
source "$project_dir"/workflow_commands/log_manager.sh #load the log manager script

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

call_function() {
	#call the cript passed in parameters, save it in logs
	#end the program if the called script has returned an error
	
	function_command=$@
	log_script_call "$container_path" "$function_command"
	$function_command #execute the command
	if [ ! $? = 0 ] #test if command failed
	then
		echo "error in $(basename $2)"
		echo "canceling dna_read"
		exit 1
	fi
}


#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

while [[ "$#" -gt 0 ]]; do
	case "$1" in
		-n_read ) if [ "$#" -lt 2 ] ; then
	    				echo "no value after param n_read" ; exit 1 ;
	    			else
	    				n_read="$2" ; shift 2 ;
	    			fi ;;
	    -h | --help ) help_function ; exit 1;;
	    -* ) echo "unknown parameter $1" ; exit 1;;
	    * ) if test -z "${container_path:-}" ; then
				container_path="${1}" ; shift ;
			elif test -z "${document_index:-}" ; then
				document_index="${1}" ; shift ;
			elif test -z "${result_path:-}" ; then
				result_path="${1}" ; shift ;
			else
				echo "unknown parameter $1" ; exit 1 ;
			fi ;;
	esac
done

if test -z "${container_path:-}"
then
	echo "container path missing"
	help_function
	exit 1
fi

if test -z "${document_index:-}"
then
	echo "document index missing"
	help_function
	exit 1
fi

if test -z "${result_path:-}"
then
	echo "result path missing"
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
    echo "error : container $container_path does not contain a document of index $document_index" 
    exit 1
fi

int_regex='^[0-9]+$'
if test "${n_read:-}" && ! [[ $n_read =~ $int_regex ]]
then
   echo "error: sequencing number ($n_read) is not a number" ; exit 1
fi


#---------------------------------------------#
######### ====== read document ====== #########
#---------------------------------------------#

is_editable=$(get_container_param $container_path "editable")

if $is_editable
then
	echo "the container is not readable"
	exit 1
fi

stored_document_path="$container_path"/$document_index


#----Molecule Selection----# 

start_time=$(date +"%s")

n_frag=$(get_doc_param $container_path $document_index "fragment_number")
start_primer=$(get_doc_param $container_path $document_index "start_primer")
stop_primer=$(get_doc_param $container_path $document_index "stop_primer")

if test -z "${n_read:-}" #if n_read has not been defined
then
	#default value of molecule sequencing number
	n_read=$(($n_frag * 50))
fi

molecule_selection_script="$project_dir"/sequencing_modules/select_molecules.py
selected_mol_path="$stored_document_path"/5_select_mol.fasta
#select molecules from container molecules with the good primers
call_function python3 "$molecule_selection_script" -i "$container_path"/container_molecules.fasta -o "$selected_mol_path" --start $start_primer --stop $stop_primer -n $n_read
mol_sel_time=$(date +"%s")


#----Sequencing & Basecalling----#

convert_fasta_script="$project_dir"/synthesis_simulation/dna_file_reader.py #script to convert fasta to fastq
simu_seq_bc_script="$project_dir"/sequencing_modules/sequencing_basecalling_simulator/sequencing_basecalling_simulator.jl
sequenced_reads_path="$stored_document_path"/6_sequenced_reads.fastq
#call_function python3 "$convert_fasta_script" "$selected_mol_path" "$sequenced_reads_path"
call_function julia "$simu_seq_bc_script" "$selected_mol_path" "$sequenced_reads_path"

seq_bc_time=$(date +"%s")


#----Consensus----#

frag_length=$(get_container_param $container_path "frag_length")
assembly_type=$(get_container_param $container_path "assembly_type")
channel_coding=$(get_doc_param $container_path $document_index "channel_coding")
consensus_path="$stored_document_path"/7_consensus.fasta

if [ $assembly_type == "indexed" ]
then
	#clustering
	clustering_script="$project_dir"/sequencing_modules/clustering/clustering #script for the clustering	
	clusters_frag_dir="$stored_document_path"/7_clusters
	
	rm -rf "$clusters_frag_dir"
	mkdir "$clusters_frag_dir"

	if $channel_coding
	then
		coded_frag_length=$(($frag_length *2)) #add redundancy from channel coding
	else
		coded_frag_length=$frag_length
	fi
		
	spacer=$(get_container_param $container_path "spacer")
	
	call_function "$clustering_script" "$sequenced_reads_path" "$clusters_frag_dir" $spacer $coded_frag_length
	
	#Consensus
	
	consensus_script="$project_dir"/sequencing_modules/cluster_consensus.py
	call_function python3 "$consensus_script" "$clusters_frag_dir" "$consensus_path" $coded_frag_length

else #ordered assembly
	start_primer=$(get_doc_param $container_path $document_index "start_primer")
	stop_primer=$(get_doc_param $container_path $document_index "stop_primer")
	
	consensus_script="$project_dir"/sequencing_modules/reads_consensus.py
	
	if $channel_coding
	then
		expected_length=$(($n_frag * $frag_length * 2)) # for a redundancy of 100%
	else
		expected_length=$(($n_frag * $frag_length)) # primers and overhangs are espected (except on last fragments ->(-4))
	fi
	
	call_function python3 "$consensus_script" -i "$sequenced_reads_path" -o "$consensus_path" --start $start_primer --stop $stop_primer -e $expected_length
fi

consensus_time=$(date +"%s")


#----Channel Decoding----#

channel_decoding_script="$project_dir"/channel_code/file_decoder.sh
decoded_sequence_path="$stored_document_path"/8_decoded_sequence.fasta

if $channel_coding
then
	# the consensus results in one unique sequence
	# it is fragmented back with the original encoded fragments length to be used in the channel decoding
	fragmented_consensus_path="$stored_document_path"/8_fragmented_consensus.fasta
	rm -rf $fragmented_consensus_path || true #remove the file if already exists
	consensus_sequence=$(head -n 2 $consensus_path | tail -1)
	frag_list=$(echo $consensus_sequence | fold -c$(($frag_length * 2)))
	for fragment in $frag_list
	do
		printf ">fragment\n$fragment\n" >> "$fragmented_consensus_path"
	done
	decoded_fragments_path="$stored_document_path"/8_decoded_fragments.fasta
	
	call_function "$channel_decoding_script" "$fragmented_consensus_path" $frag_length "$decoded_fragments_path" "$container_path"/$document_index/validity_check.txt
	
	#the decoded fragments are reassembled in an unique sequence
	#select even lines and delete all \n 
	printf ">decoded_consensus\n" > "$decoded_sequence_path"
	sed -n '0~2p' $decoded_fragments_path | tr --delete '\n' >> "$decoded_sequence_path"
else
	call_function cp "$consensus_path" "$decoded_sequence_path"
fi

channel_time=$(date +"%s")


#----Source Decoding----#

encoding_data=$(get_all_doc_param $container_path $document_index)

source_decoding_script="$project_dir"/source_encoding/source_decoding.py
call_function python3 "$source_decoding_script" -i "$decoded_sequence_path" -o "$result_path" -t "$assembly_type" -l $frag_length --data "$encoding_data"

echo "Document $document_index of $container_path successfully saved to $result_path !"

end_time=$(date +"%s")
times_file="$stored_document_path"/workflow_times.txt
echo "dna_read : $(($end_time - $start_time)) s" >> "$times_file"
echo "   > molecules_select       : $(($mol_sel_time - $start_time)) s" >> "$times_file"
echo "   > sequencing_basecalling : $(($seq_bc_time - $mol_sel_time)) s" >> "$times_file"
echo "   > consensus              : $(($consensus_time - $seq_bc_time)) s" >> "$times_file"
echo "   > channel_decoding       : $(($channel_time - $consensus_time)) s" >> "$times_file"
echo "   > source_decoding        : $(($end_time - $channel_time)) s" >> "$times_file"

exit 0
