#!/bin/bash
set -u #exit and display error message if a variable is empty 

#DNA WORFLOW - run the complete cycle of archiving and extracting a document in a container
# dna_workflow [-no_read] Dname Cname", or also dna_workflow [-no_read] for default container and document

project_dir="$(dirname $0)/.." #parent of the directory containing this script
commands_dir="$project_dir"/workflow_commands
source "$project_dir"/workflow_commands/log_manager.sh #load the log manager script

call_function() {
	#call the cript passed in parameters, save it in logs
	#end the program if the called script has returned an error
	
	function_command=$@
	log_script_call "$container_path"/log_file.log "$function_command"
	$function_command #execute the command
	if [ ! $? = 0 ] #test if command failed
	then
		echo "error in $(basename $1)"
		echo "canceling dna_workflow"
		exit 1
	fi
}

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

case "$1" in
    -no_read ) no_read=true ; document_path="${2}"; container_path="${3}";;
    -* ) echo "unknown parameter $1" ; exit 1;;
    * ) no_read=false ; document_path="${1}"; container_path="${2}" ;;
esac

if test -z "$document_path"
then
	document_path="documents_test/doc.txt" #img_200.png"
fi

if test -z "$container_path"
then
	container_path="test_workflow"
fi

#------------------------------------------------#
######### ====== run the workflow ====== #########
#------------------------------------------------#

if [ -d "$container_path" ]
then
	rm -rf "$container_path"_old
	mv "$container_path" "$container_path"_old #save the previous workflow container
fi

call_function "$commands_dir"/dna_create.sh -sim -nocd -fl 60 "$container_path" #TODO fl

call_function "$commands_dir"/dna_add.sh "$document_path" "$container_path"

call_function "$commands_dir"/dna_store.sh "$container_path"

if $no_read #skip the reading part
then
	exit 0
fi

document_name="$(basename $document_path)"
call_function "$commands_dir"/dna_read.sh "$container_path" 0 "$container_path"/0/result_"$document_name"


#command to copy the total workflow from the dnarxiv server
# scp -r oboulle@dnarxiv.irisa.fr:~/Documents/workflow_commands/test_workflow .
