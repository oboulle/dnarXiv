#!/bin/bash

#nohup ./multiple_call.sh &

SCRIPT_PATH="workflow_4.sh"

for i in 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000
do
	./workflow_4.sh -n "$i ";
done


#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
