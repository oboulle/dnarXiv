#!/bin/bash
set -u #exit and display error message if a variable is empty 

#DNA ADD - add a document to the container

project_dir="$(dirname $(dirname ${BASH_SOURCE}))" #parent of the directory containing this script
source "$project_dir"/workflow_commands/metadata_manager.sh #load the xml manager script
source "$project_dir"/workflow_commands/log_manager.sh #load the log manager script

help_function() {
   echo ""
   echo "Usage: dna_add [-nocd] CName DName "
   echo -e "\t-nocd : turn off channel encoding"
   echo -e "\tCName : path to the container"
   echo -e "\tDName : path to the document"
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
		echo "canceling dna_add"
		rm -rf $stored_document_path
		del_document "$container_path" $doc_index
		exit 1
	fi
}

#----- default parameters -----#
channel_coding=true

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

while [[ "$#" -gt 0 ]]; do
	case "$1" in
		-nocd ) channel_coding=false ; shift ;;
		-h | --help ) help_function ; exit 1 ;;
		-* ) echo "unknown parameter $1" ; exit 1 ;;
		* ) if test -z "${container_path:-}" ; then
				container_path="${1}" ; shift ;
			elif test -z "${document_path:-}" ; then
				document_path="${1}" ; shift ;
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

if test -z "${document_path:-}"
then
	echo "document path missing"
	help_function
	exit 1
fi

if [ ! -d "$container_path" ] 
then
    echo "error : container $container_path not found" 
    exit 1
fi

if [ ! -f "$document_path" ] 
then
    echo "error : document $document_path not found" 
    exit 1
fi 

#--------------------------------------------#
######### ====== add document ====== #########
#--------------------------------------------#

time=$(date +"%s")

is_editable=$(get_container_param $container_path "editable")

if ! $is_editable
then
	echo "the container is not editable"
	exit 1
fi

# add a new document in the xml file and get its new index
doc_index=$(add_document "$container_path")

stored_document_path="$container_path"/$doc_index
 
mkdir -p "$stored_document_path"

#----Source Encoding----#

# get the fragment length from the container options
frag_length=$(get_container_param $container_path "frag_length")

assembly_type=$(get_container_param $container_path "assembly_type")

source_encoding_script="$project_dir"/source_encoding/source_encoding.py
source_path="$stored_document_path"/0_encoded_source.fasta

call_function python3 "$source_encoding_script" -i "$document_path" -o "$source_path" -l $frag_length -t "$assembly_type" --cont_doc $stored_document_path


#----Channel Encoding----#
#TODO channel coding for whole sequence
channel_encoding_script="$project_dir"/channel_code/encode_from_file.jl 
channel_path="$stored_document_path"/1_channel.fasta

add_doc_param "$container_path" $doc_index "channel_coding" $channel_coding

if $channel_coding
then
	call_function julia "$channel_encoding_script" "$source_path" "$channel_path"
else
	call_function cp "$source_path" "$channel_path"
fi



echo "Document $(basename $document_path) successfully added to container $(basename $container_path) !"

times_file="$stored_document_path"/workflow_times.txt

echo workflow $(date +"%Hh%Mm%S") >> "$times_file"

end_time=$(date +"%s")
echo "dna_add : $(($end_time - $time)) s" >> "$times_file"
exit 0
