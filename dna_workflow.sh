#!/bin/bash

#DNA WORFLOW - run the complete cycle of archiving and extracting a document in a container

check_error_function () { #end the program if the previously called script has returned an error
	if [ ! $? = 0 ]
	then
		echo "error in $1"
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
document_name="img_50.png"

rm -rf "$container_name"_old
mv "$container_name" "$container_name"_old #save the previous workflow

"$commands_dir"/dna_create.sh -sim -fl 100 "$container_name"
check_error_function "dna_create"

"$commands_dir"/dna_add.sh "$commands_dir"/"$document_name" "$container_name"
check_error_function "dna_add"

"$commands_dir"/dna_store.sh "$container_name"
check_error_function "dna_store"

"$commands_dir"/dna_read.sh "$container_name" 0 "$container_name"/result_"$document_name"
check_error_function "dna_read"


#----Result Analysis----#

save_results_script="$commands_dir"/result_analysis/save_results.py
source_path="$container_name"/0/0_source.fasta
reconstructed_source_path="$container_name"/0/11_reconstructed_source.fasta
n_seq=$(($(wc -l < "$container_name"/0/6_select_mol.fasta)/2))
times_file="$container_path"/workflow_times.txt
results_file="$container_name"/results.txt
python3 "$save_results_script" "$source_path" "$reconstructed_source_path" $n_seq "$times_file" "$results_file"


exit 0

#command to copy the total workflow from the dnarxiv server
# scp -r oboulle@dnarxiv.irisa.fr:~/Documents/workflow_commands/test_workflow .
