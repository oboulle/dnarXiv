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

commands_dir="$project_dir"/workflow_commands
container_name="test_workflow"
rm -rf "$container_name"
"$commands_dir"/dna_create.sh -sim -fl 100 "$container_name" 
"$commands_dir"/dna_add.sh "$commands_dir"/doc.txt "$container_name"
"$commands_dir"/dna_store.sh "$container_name" 
"$commands_dir"/dna_list.sh "$container_name"
"$commands_dir"/dna_read.sh "$container_name" 0 resultat_du_workflow.txt

exit 0