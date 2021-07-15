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

#get metadata informations of META
#if META not already sequenced
#sequence META dir

#get metadata informations of DI

#sequence DI

echo "Document $document_index of $container_path successfully saved to $document_path !"
exit 0
