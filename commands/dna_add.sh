#!/bin/bash

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

help_function() {
   echo ""
   echo "Usage: dna_add Dname Cname"
   exit 1 # Exit script after printing help
}

document_path="${1}"
container_path="${2}"

if [ "$#" != 2 ]
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

if [ ! -f "$document_path" ] 
then
    echo "error : document $document_path not found" 
    exit 1
fi 

#--------------------------------------------#
######### ====== add document ====== #########
#--------------------------------------------#

cdi_file="$container_path/.cdi"
document_index=$(head -n 1 "$cdi_file")

if (( $document_index < 0 ))
then
	echo "the container is not editable"
	exit 1
fi

stored_document_path="$container_path/$document_index"
 
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
document_index $document_index
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

channel_encoding_script="$project_dir/synthesis_simulation/channel_encoding/channel_encoding.py" 
channel_path="$stored_document_path/channel.fasta"

python3 $channel_encoding_script "$source_path" "$channel_path" #TODO
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
	echo "error in homopolymere deletion"
	#cancel_dna_add#TODO
	#exit 1
fi

#----Update .cdi----#
document_index=$((document_index+1))

cat > $cdi_file << eof
$document_index
eof


echo "Document $document_path successfully added to container $container_path !"
exit 0
