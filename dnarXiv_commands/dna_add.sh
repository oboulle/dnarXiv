#!/bin/bash

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

help_function() {
   echo ""
   echo "Usage: dna_add Dname Cname"
   echo -e "\tDname : path to the document"
   echo -e "\tCname : path to the container"
   echo ""
   exit 1 # Exit script after printing help
}

while true; do
  case "$1" in
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
cdi_file="$container_path/.cdi"
container_index=$(head -n 1 "$cdi_file")

if (( $container_index < 0 ))
then
	echo "the container is not editable"
	exit 1
fi

stored_document_path="$container_path/$container_index"
 
mkdir -p "$stored_document_path"

if [ $HOME == "/Users/oboulle" ]
then
  project_dir="/Users/oboulle/Documents"
else
  project_dir="/home/oboulle/Documents"
fi

meta_file="$stored_document_path/.meta"
cat > $meta_file << eof
document_index $container_index
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

source_encoding_script="$project_dir/source_encoding/source_encoding.py" 
source_path="$stored_document_path/source.fasta"

python3 $source_encoding_script "$document_path" "$source_path" "$frag_length" "$meta_file" #TODO
if [ ! $? = 0 ]
then
	echo "error in source encoding"
	cancel_dna_add#TODO
	exit 1
fi

#----Channel Encoding----#

channel_encoding_script="$project_dir/channel_code/encode_from_file.jl" 
channel_path="$stored_document_path/channel.fasta"

julia $channel_encoding_script "$source_path" "$channel_path" #TODO
if [ ! $? = 0 ]
then
	echo "error in channel encoding"
	cancel_dna_add#TODO
	exit 1
fi

#----Homopolymere Deletion----#

h_deletion_script="$project_dir/synthesis_simulation/homoplymere_deletion/homopolymere_deletion.py" 
fragments_path="$stored_document_path/fragments.fasta"

echo "homopolymere deletion not implemented yet; skipping..."
cp $channel_path $fragments_path

: 'python3 $h_deletion_script "$channel_path" "$fragments_path" #TODO
if [ ! $? = 0 ]
then
	#echo "error in homopolymere deletion"
	#cancel_dna_add#TODO
	#exit 1
fi
'
#----Update .cdi----#

container_index=$((container_index+1))

cat > $cdi_file << eof
$container_index
eof


echo "Document $document_path successfully added to container $container_path !"
end_time=$(date +"%s")
echo "dna_add : $(($end_time - $time)) s" >> workflow_times.txt
exit 0
