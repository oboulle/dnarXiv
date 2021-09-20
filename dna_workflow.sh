#!/bin/bash

#DNA WORFLOW - run the complete cycle of archiving and extracting a document in a container
# dna_workflow [-no_read] Cname Dname", or also dna_workflow [-no_read] for default container and document

check_error_function () { #end the program if the previously called script has returned an error
	if [ ! $? = 0 ]
	then
		echo "error in $1"
		save_results_function
		exit 1
	fi
}

if [ $HOME == "/Users/oboulle" ]
then
	project_dir="/Users/oboulle/Documents"
elif [ $HOME == "/udd/oboulle" ]
then
	project_dir="/udd/oboulle/Documents"
else
	project_dir="/home/oboulle/Documents"
fi

commands_dir="$project_dir"/workflow_commands

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

case "$1" in
    -no_read ) no_read=true ; document_name="${2}"; container_name="${3}";;
    -* ) echo "unknown parameter $1" ; exit 1;;
    * ) no_read=false ; document_name="${1}"; container_name="${2}" ;;
esac

if test -z "$document_name"
then
	document_name="doc.txt" #"img_50.png"
fi

if test -z "$container_name"
then
	container_name="test_workflow"
fi

#------------------------------------------------#
######### ====== run the workflow ====== #########
#------------------------------------------------#

rm -rf "$container_name"_old
mv "$container_name" "$container_name"_old #save the previous workflow

"$commands_dir"/dna_create.sh -sim "$container_name"
check_error_function "dna_create"

"$commands_dir"/dna_add.sh "$commands_dir"/"$document_name" "$container_name"
check_error_function "dna_add"

"$commands_dir"/dna_store.sh "$container_name"
check_error_function "dna_store"

if $no_read #skip the reading part
then
	exit 0
fi

"$commands_dir"/dna_read.sh "$container_name" 0 "$container_name"/result_"$document_name"
check_error_function "dna_read"


#command to copy the total workflow from the dnarxiv server
# scp -r oboulle@dnarxiv.irisa.fr:~/Documents/workflow_commands/test_workflow .
