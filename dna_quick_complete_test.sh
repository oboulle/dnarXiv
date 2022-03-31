#!/bin/bash

if [ $HOME == "/Users/oboulle" ]
then
	project_dir="/Users/oboulle/Documents"
elif [ $HOME == "/udd/oboulle" ]
then
	project_dir="/udd/oboulle/Documents"
else
	project_dir="/home/oboulle/Documents"
fi

check_error_function () { 
	if [ ! $? = 0 ]
	then
		echo "error in $1"
		echo "cancel dna_add"
		rm -rf $stored_document_path
		exit 1
	fi
}

commands_dir="$project_dir"/workflow_commands
container_name="test_workflow"
rm -rf "$container_name"

"$commands_dir"/dna_create.sh -sim -fl 100 "$container_name" 
check_error_function "dna_create"

"$commands_dir"/dna_add.sh -nocd "$commands_dir"/documents_test/doc.txt "$container_name"
check_error_function "dna_add"

"$commands_dir"/dna_add.sh "$commands_dir"/documents_test/img_50.png "$container_name"
check_error_function "dna_add"

"$commands_dir"/dna_store.sh "$container_name"
check_error_function "dna_store"

"$commands_dir"/dna_list.sh "$container_name"
check_error_function "dna_list"

"$commands_dir"/dna_read.sh "$container_name" 0 resultat_du_workflow.txt
check_error_function "dna_read"


exit 0