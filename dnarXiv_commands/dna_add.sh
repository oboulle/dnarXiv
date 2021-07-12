#!/bin/bash

is_meta=false
#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

help_function() {
   echo ""
   echo "Usage: dna_add [-meta] Dname Cname"
   echo -e "\t-meta : used to add a metadata file to the container"
   echo -e "\tDname : path to the document"
   echo -e "\tCname : path to the container"
   echo ""
   exit 1 # Exit script after printing help
}

while true; do
  case "$1" in
    -meta ) is_meta=true ; shift ;;
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

cdi_file="$container_path/.cdi"
container_index=$(head -n 1 "$cdi_file")

if (( $container_index < 0 ))
then
	echo "the container is not editable"
	exit 1
fi

if [ $is_meta = true ]
then
	container_index="META"
fi

stored_document_path="$container_path/$container_index"
 
mkdir -p "$stored_document_path"

working_dir=$(pwd)

if [ $HOME == "/Users/oboulle" ]
then
  project_dir="/Users/oboulle/Documents"
else
  project_dir="/home/oboulle/Documents"
fi

meta_file="$stored_document_path/.meta"
cat > $meta_file << eof
container_index $container_index
creation_date $(date +'%d/%m/%Y %R')
eof

cancel_dna_add() {
	echo "cancel dna_add"
	rm -rf $stored_document_path
}

#----Source Encoding----#

# get the fragment length from the container options
while read var value; do
    export "$var"="$value"
done < "$container_path/.options"

source_encoding_script="$project_dir/synthesis_simulation/source_encoding/source_encoding.py" 
source_path="$stored_document_path/source.fasta"

python3 $source_encoding_script "$document_path" "$source_path" "$frag_length" "$meta_file" #TODO
if [ ! $? = 0 ]
then
	echo "error in source encoding"
	#cancel_dna_add#TODO
	#exit 1
fi

#----Channel Encoding----#

channel_encoding_script="$project_dir/channel_coding/LDPC/encode_from_file.jl" 
channel_path="$stored_document_path/channel.fasta"

julia $channel_encoding_script "$source_path" "$channel_path" #TODO
if [ ! $? = 0 ]
then
	echo "error in channel encoding"
	#cancel_dna_add#TODO
	#exit 1
fi

#----Homopolymere Deletion----#

h_deletion_script="$project_dir/synthesis_simulation/homoplymere_deletion/homopolymere_deletion.py" 
fragments_path="$stored_document_path/fragments.fasta"

python3 $h_deletion_script "$channel_path" "$fragments_path" #TODO
if [ ! $? = 0 ]
then
	#echo "error in homopolymere deletion"
	echo "homopolymere deletion not implemented yet; skipping..."
	#cancel_dna_add#TODO
	#exit 1
fi

#----Update .cdi----#

#no update if the document is the meta file of the container
if [ $is_meta = true ]
then
	echo "Meta document successfully added to container $container_path !"
	exit 0
fi
container_index=$((container_index+1))

cat > $cdi_file << eof
$container_index
eof


echo "Document $document_path successfully added to container $container_path !"
exit 0
