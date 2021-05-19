#!/bin/bash

#nohup ./multiple_call.sh &

SCRIPT_PATH="workflow_4.sh"
n_iter=2

for (( i=0; i <= $n_iter; i++ ));
do 
	echo "_________run $i _________";
	. "$SCRIPT_PATH" -n 125
	. "$SCRIPT_PATH" -n 150
	. "$SCRIPT_PATH" -n 175
	. "$SCRIPT_PATH" -n 200
	. "$SCRIPT_PATH" -n 225
	. "$SCRIPT_PATH" -n 250
	. "$SCRIPT_PATH" -n 275
	. "$SCRIPT_PATH" -n 300
	. "$SCRIPT_PATH" -n 325
	. "$SCRIPT_PATH" -n 350
	. "$SCRIPT_PATH" -n 375
	. "$SCRIPT_PATH" -n 400
done



#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
