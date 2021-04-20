#!/bin/bash


SCRIPT_PATH="workflow_4.sh"

n_iter=10

for (( i=0; i <= $n_iter; i++ ));
do 
	echo "_________run $i _________";
	. "$SCRIPT_PATH"; 
done



#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
