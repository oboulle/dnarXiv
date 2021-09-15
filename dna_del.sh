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
document_index="${2}"

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

if [ ! -d "$container_path"/$document_index ] 
then
    echo "error : $container_path does not contain a document of index $document_index" 
    exit 1
fi

#-----------------------------------------------#
######### ====== delete document ====== #########
#-----------------------------------------------#

cdi_file="$container_path"/.cdi
current_document_index=$(head -n 1 "$cdi_file")

if (( $current_document_index < 0 ))
then
	echo "the container is not editable (the documents have been stored)"
	exit 1
fi

rm -rf "$container_path"/$document_index
echo "Document $document_index successfully deleted from $container_path !"

exit 0
