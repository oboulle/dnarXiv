#!/bin/bash


n_mol_list=`seq 5 5 40` #'seq first step last'
same_script_count=3 #number of time the script is executed with the same parameters
container_name=test_perf_container
document_name=doc.txt

(exec ../dna_workflow.sh -no_read $document_name $container_name) #init a container 

for n_mol in $n_mol_list; do
	for i in $(seq $same_script_count); do
		echo ______________ n_mol $n_mol : $i ______________
		rm -rf $container_name/0/{6*,7*,8*,9*,10*,11*,12*,workflow_times.txt} #remove previous reading files
	    (exec ../dna_read.sh -n_mol $n_mol $container_name 0 results/${n_mol}_${i}_$document_name)
	    python3 save_results.py $container_name/0 results_workflow.txt
	done
done

#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
