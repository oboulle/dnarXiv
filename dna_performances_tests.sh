#!/bin/bash


SCRIPT_PATH="./multiple_call.sh"

(cd "workflow_10k" && exec "$SCRIPT_PATH")
(cd "workflow_25k" && exec "$SCRIPT_PATH")
(cd "workflow_50k" && exec "$SCRIPT_PATH")
(cd "workflow_100k" && exec "$SCRIPT_PATH")

#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
