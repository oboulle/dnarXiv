#!/bin/bash

#DNA LIST - display the informations of every documents in the container

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

source ./metadata_manager.sh #load the xml manager script
meta_file="$container_path"/metadata.xml

if [[ $(get_container_param $meta_file "number_of_documents") -eq 0 ]]
then
	echo "this container is empty"
	exit 0
fi

doc_params_list=(creation_date doc_name doc_type color height width fragment_number channel_coding starting_sequence start_primer stop_primer)

for dir_name in $(find "$container_path"/* -maxdepth 0 -type d -printf "%f\n") ; do
	#read the metadata for each stored document and asign it to variables
	printf "$dir_name :\n"
	for doc_param in "${doc_params_list[@]}" ; do
		param_value=$(get_doc_param $meta_file $dir_name $doc_param)
		if test "$param_value"
		then
			printf "\t$doc_param $param_value\n"
		fi
	done
done

exit 0
