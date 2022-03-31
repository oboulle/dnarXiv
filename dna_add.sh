#!/bin/bash

#DNA ADD - add a document to the container

help_function() {
   echo ""
   echo "Usage: dna_add Dname Cname [-nocd]"
   echo -e "\tDname : path to the document"
   echo -e "\tCname : path to the container"
   echo -e "\t-nocd : turn off channel encoding"
   echo ""
   exit 1 # Exit script after printing help
}

# end the program if the previously called script has returned an error
check_error_function () { 
	if [ ! $? = 0 ]
	then
		echo "error in $1"
		echo "cancel dna_add"
		rm -rf $stored_document_path
		del_document "$meta_file" $doc_index
		exit 1
	fi
}

#----- default parameters -----#
channel_coding=true

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#


while true; do
  case "$1" in
	-nocd ) channel_coding=false ; shift ;;
    -h | --help ) help_function ; exit 1;;
    -* ) echo "unknown parameter $1" ; exit 1;;
    * ) document_path="${1}" ; container_path="${2}" ; break ;;
  esac
done

if test -z "$container_path"
then
	echo "container path missing"
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

source ./metadata_manager.sh #load the xml manager script
meta_file="$container_path"/metadata.xml


is_editable=$(get_container_param $meta_file "editable")

if ! $is_editable
then
	echo "the container is not editable"
	exit 1
fi

#get list of current documents in the container
doc_list=$(find "$container_path"/* -maxdepth 0 -type d -printf "%f\n")

#if none, new document of index 0
if [ -z "$doc_list" ]; then
	doc_index=0
else
#else index equals to max index +1
	max_index=$(echo "${doc_list[*]}" | sort -nr | head -n1)
	doc_index=$((max_index+1))
fi

stored_document_path="$container_path"/$doc_index
 
mkdir -p "$stored_document_path"

if [ $HOME == "/Users/oboulle" ] #TODO exportable
then
	project_dir="/Users/oboulle/Documents"
elif [ $HOME == "/udd/oboulle" ]
then
	project_dir="/udd/oboulle/Documents"
else
	project_dir="/home/oboulle/Documents"
fi

add_document "$meta_file" $doc_index

# get the fragment length from the container options
frag_length=$(get_container_param $meta_file "frag_length")

#----Source Encoding----#

source_encoding_script="$project_dir"/source_encoding/source_encoding.py
source_path="$stored_document_path"/1_fragments.fasta

python3 "$source_encoding_script" -i "$document_path" -o "$source_path" -l $frag_length --meta "$meta_file" --doc_index $doc_index
check_error_function "source encoding"


#----Channel Encoding----#

channel_encoding_script="$project_dir"/channel_code/encode_from_file.jl 
channel_path="$stored_document_path"/2_channel.fasta

if $channel_coding
then
	"$channel_encoding_script" "$source_path" "$channel_path"
else
	cp "$source_path" "$channel_path"
fi

add_doc_param "$meta_file" $doc_index "channel_coding" $channel_coding
check_error_function "channel encoding"


#----Homopolymere Deletion----#

#h_deletion_script="$project_dir"/synthesis_simulation/homoplymere_deletion/homopolymere_deletion.py
fragments_path="$stored_document_path"/3_final_fragments.fasta

echo "homopolymere deletion not implemented yet; skipping..."
cp "$channel_path" "$fragments_path"

: 'python3 $h_deletion_script "$channel_path" "$fragments_path" #TODO
check_error_function "homopolymere deletion"
'

echo "Document $document_path successfully added to container $container_path !"

times_file="$stored_document_path"/workflow_times.txt

echo workflow $(date +"%Hh%Mm%S") >> "$times_file"

end_time=$(date +"%s")
echo "dna_add : $(($end_time - $time)) s" >> "$times_file"
exit 0
