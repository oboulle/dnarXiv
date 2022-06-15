#!/bin/bash
set -u #exit and display error message if a variable is empty 

#DNA LIST - display the informations of every documents in the container

project_dir="$(dirname $(dirname ${BASH_SOURCE}))" #parent of the directory containing this script
source "$project_dir"/workflow_commands/metadata_manager.sh #load the xml manager script

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

help_function() {
   echo ""
   echo "Usage: dna_list CName"
   echo -e "\tCName : path to the container"
   echo ""
   exit 1 # Exit script after printing help
}


if [ "$#" != 1 ]
then
	echo "error invalid number of arguments"
	help_function
	exit 1
else
	container_path="${1}"
fi

if [ ! -d "$container_path" ] 
then
    echo "error : container $container_path not found" 
    exit 1
fi

#----------------------------------------------#
######### ====== list documents ====== #########
#----------------------------------------------#

if [[ $(get_container_param $container_path "number_of_documents") -eq 0 ]]
then
	echo "this container is empty"
	exit 0
fi

doc_params_list=(creation_date doc_name doc_type color height width fragment_number channel_coding start_sequence start_primer stop_primer)

#loop over names of directories in the container (should only be "0","1","2",...)
for dir_name in $(find "$container_path" -maxdepth 1 -mindepth 1 -type d -exec basename {} ';') ; do
	#read the metadata for each stored document and asign it to variables
	printf "$dir_name :\n"
	for doc_param in "${doc_params_list[@]}" ; do
		param_value=$(get_doc_param $container_path $dir_name $doc_param)
		if [ "$param_value" != "undefined" ]
		then
			printf "\t$doc_param $param_value\n"
		fi
	done
done

exit 0
