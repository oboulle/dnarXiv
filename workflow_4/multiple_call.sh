#!/bin/bash

#nohup ./multiple_call.sh &

SCRIPT_PATH="workflow_4.sh"
SCRIPT_PATH_2="workflow_4_C.sh"
n_iter=10

for (( i=0; i <= $n_iter; i++ ));
do 
	echo "_________run $i _________";
	. "$SCRIPT_PATH";
	. "$SCRIPT_PATH_2";
done



#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
