#!/bin/bash

#DNA WORFLOW - run the complete cycle of archiving and extracting a document in a container

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
container_name="test_workflow"
document_name="doc.txt" #"img_50.png"

save_results_function () {
	save_results_script="$commands_dir"/result_analysis/save_results.py
	results_file="$commands_dir"/results_workflow.txt
	python3 "$save_results_script" "$container_name"/0 "$results_file"
}

rm -rf "$container_name"_old
mv "$container_name" "$container_name"_old #save the previous workflow

"$commands_dir"/dna_create.sh -sim "$container_name"
check_error_function "dna_create"

"$commands_dir"/dna_add.sh "$commands_dir"/"$document_name" "$container_name"
check_error_function "dna_add"

"$commands_dir"/dna_store.sh "$container_name"
check_error_function "dna_store"

if test "$1" #number of molecules to read is defined
then
	"$commands_dir"/dna_read.sh -n_mol $1 "$container_name" 0 "$container_name"/result_"$document_name"
else
	"$commands_dir"/dna_read.sh "$container_name" 0 "$container_name"/result_"$document_name"
fi

check_error_function "dna_read"

save_results_function


#command to copy the total workflow from the dnarxiv server
# scp -r oboulle@dnarxiv.irisa.fr:~/Documents/workflow_commands/test_workflow .
