#!/bin/bash


n_mol_list=`seq 5 5 40` #'seq first step last'
container_name=test_perf_container
document_name=doc.txt

(exec ../dna_workflow.sh -no_read $document_name $container_name) #init a container 

for n_mol in $n_mol_list; do
	echo ______________ n_mol $n_mol : $i ______________
	rm -rf $container_name/0/{6*,7*,8*,9*,10*,11*,12*} #remove previous reading files
	mv $container_name/0/workflow_times.txt $container_name/workflow_times.txt
	(exec ../dna_read.sh -n_mol $n_mol $container_name 0 results/${n_mol}_${i}_$document_name)
	python3 save_results.py $container_name/0 results_workflow.txt
done

#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
