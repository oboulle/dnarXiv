#!/bin/bash

#DNA DEL - remove a document from the container

help_function() {
   echo ""
   echo "Usage: dna_del Cname DI"
   echo -e "\tCname : path to the container"
   echo -e "\tDI : index of the document to delete"
   echo ""
   exit 1 # Exit script after printing help
}

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

container_path="${1}"
doc_index="${2}"

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

if [ ! -d "$container_path"/$doc_index ] 
then
    echo "error : $container_path does not contain a document of index $doc_index" 
    exit 1
fi

#-----------------------------------------------#
######### ====== delete document ====== #########
#-----------------------------------------------#

project_dir="$(dirname $0)/.." #parent of the directory containing this script

source "$project_dir"/workflow_commands/metadata_manager.sh #load the xml manager script
meta_file="$container_path"/metadata.xml

is_editable=$(get_container_param $meta_file "editable")

if ! $is_editable
then
	echo "the container is not editable"
	exit 1
fi

rm -rf "$container_path"/$doc_index
del_document "$meta_file" $doc_index

echo "Document $doc_index successfully deleted from $container_path !"

exit 0
