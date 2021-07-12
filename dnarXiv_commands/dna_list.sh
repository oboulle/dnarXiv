#!/bin/bash

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

help_function() {
   echo ""
   echo "Usage: dna_list Cname"
   echo -e "\tCname : path to the container"
   echo ""
   exit 1 # Exit script after printing help
}

container_path="${1}"

if [ "$#" != 1 ]
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

#--------------------------------------------#
######### ====== list documents ====== #########
#--------------------------------------------#

if [[ $(find $container_path -maxdepth 1 -type d | wc -l) -eq 1 ]]
then
	echo "this container is empty"
	exit 0
fi

for directory in $container_path/*/ ; do
	#read the data in the .meta file of each stored document and asign it to variables
    while read var value; do
    	export "$var"="$value"
	done < "$directory/.meta"
	printf "$document_index :\n"
	printf "\tcreation date $creation_date\n"
	printf "\tname $doc_name\n"
	printf "\ttype $type\n"
	printf "\t$n_frag fragments\n"
done

exit 0
