#!/bin/bash

if [ $HOME == "/Users/oboulle" ]
then
  project_dir="/Users/oboulle/Documents"
else
  project_dir="/home/oboulle/Documents"
fi

commands_dir="$project_dir/workflow_global/dnarXiv_commands"
rm -rf test_workflow
$commands_dir/dna_create.sh -sim -fl 100 test_workflow 
$commands_dir/dna_add.sh $commands_dir/doc.txt test_workflow
$commands_dir/dna_store.sh test_workflow 
$commands_dir/dna_list.sh test_workflow
$commands_dir/dna_read.sh test_workflow 0 resultat_du_workflow.txt

exit 0