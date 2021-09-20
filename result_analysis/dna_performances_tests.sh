#!/bin/bash

SCRIPT_PATH="../dna_workflow.sh"
n_mol_list=`seq 40 -5 0` #'seq first step last'; decrement because lowest n_mol will be more likely to fail
same_script_count=5 #number of time the exact same script is executed

for n_mol in $n_mol_list; do
	for i in $(seq $same_script_count); do
		echo ______________ n_mol $n_mol : $i ______________
	    (exec "$SCRIPT_PATH" $n_mol)
	done
done

#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
