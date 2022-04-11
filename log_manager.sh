#!/bin/bash
set -u #exit and display error message if a variable is empty 

#script for log management in the containers

log_file_name="log_file.log" #name of the log file in containers

log_script_call() {
	#register the call of a script
	#1 container path
	#2 call of the script ex: "./script.py param1 param2"
	log_file="$1"/$log_file_name
	script_call="$2"
	
	#remove the path to the script name
	script_name=$(basename $(echo "$script_call" | head -n1 | cut -d " " -f1))
	script_params=$(echo "$script_call" | cut -d ' ' -f 2-)
	
	current_date=$(date '+%Y-%m-%d %H:%M:%S')
	
	echo "[$current_date] [SCRIPT] [$script_name $script_params]" >> "$log_file"
}

log_metadata_call() {
	#register the call to the metadata manager
	#1 container path
	#2 call of the metadata manager ex: "change_param param1 param2"
	log_file="$1"/$log_file_name
	meta_call="${@:2}"
	
	current_date=$(date '+%Y-%m-%d %H:%M:%S')
	
	echo "[$current_date] [METADATA] [$meta_call]" >> "$log_file"
}
