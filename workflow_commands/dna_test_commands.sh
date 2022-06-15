#!/bin/bash
set -u #exit and display error message if a variable is empty 

project_dir="$(dirname ${BASH_SOURCE})/.." #parent of the directory containing this script
commands_dir="$project_dir"/workflow_commands
source "$project_dir"/workflow_commands/log_manager.sh #load the log manager script


call_function() {
	#call the cript passed in parameters, save it in logs
	#end the program if the called script has returned an error
	
	function_command=$@
	log_script_call "$container_path" "$function_command"
	$function_command #execute the command
	if [ ! $? = 0 ] #test if command failed
	then
		echo "error in $(basename $1)"
		echo "canceling quick complete test"
		exit 1
	fi
}


container_path="test_workflow"
rm -rf "$container_path"

"$commands_dir"/dna_create.sh -sim -fl 60 "$container_path" 

call_function "$commands_dir"/dna_add.sh -nocd "$container_path" "$commands_dir"/documents_test/doc_10_frag.txt #dij.txt

#call_function "$commands_dir"/dna_add.sh -nocd "$container_path" "$commands_dir"/documents_test/img_50.png

call_function "$commands_dir"/dna_store.sh "$container_path"

call_function "$commands_dir"/dna_list.sh "$container_path"

call_function "$commands_dir"/dna_read.sh "$container_path" 0 $container_path/resultat_0.txt

#call_function "$commands_dir"/dna_read.sh "$container_path" 1 $container_path/resultat_1.png

#command to copy the total workflow from the dnarxiv server
# scp -r oboulle@dnarxiv.irisa.fr:~/Documents/workflow_commands/test_workflow .

exit 0