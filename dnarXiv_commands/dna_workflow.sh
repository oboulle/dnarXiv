#!/bin/bash

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

commands_dir="$project_dir/workflow_global/dnarXiv_commands"
rm -rf test_workflow_old
mv test_workflow test_workflow_old
$commands_dir/dna_create.sh -sim -fl 100 test_workflow
check_error_function "dna_create"
$commands_dir/dna_add.sh $commands_dir/doc.txt #test_workflow
check_error_function "dna_add"
#$commands_dir/dna_add.sh $commands_dir/img_small.png test_workflow
$commands_dir/dna_store.sh test_workflow 
check_error_function "dna_store"
#$commands_dir/dna_list.sh test_workflow
$commands_dir/dna_read.sh test_workflow 0 test_workflow/resultat_du_workflow.txt
check_error_function "dna_read"

exit 0

# scp -r oboulle@dnarxiv.irisa.fr:~/Documents/workflow_global/dnarXiv_commands/test_workflow/0 .