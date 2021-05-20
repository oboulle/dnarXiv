#!/bin/bash

#nohup ./multiple_call.sh &

SCRIPT_PATH="workflow_4.sh"

for i in 125 150 175 200 225 250 275 300 325 350 375 400
do
	./workflow_4.sh -n "$i ";
done


#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
