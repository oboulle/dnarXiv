#!/bin/bash
set -u #exit and display error message if a variable is empty 

#DNA STORE - synthetise the documents of the container into molecules

project_dir="$(dirname $(dirname ${BASH_SOURCE}))" #parent of the directory containing this script
source "$project_dir"/workflow_commands/metadata_manager.sh #load the xml manager script
source "$project_dir"/workflow_commands/log_manager.sh #load the log manager script

help_function() {
   echo ""
   echo "Usage: dna_store [-n_synth int] [-i_error float] [-d_error float] [-s_error float] CName"
   echo -e "\t-n_synth : specify number of synthesis per fragments [default = 100]"
   echo -e "\t-i_error : specify insertion error rate in synthesis simulation [default = 0]"
   echo -e "\t-d_error : specify deletion error rate in synthesis simulation [default = 0]"
   echo -e "\t-s_error : specify substitution error rate in synthesis simulation [default = 0]"
   echo -e "\tCName : path to the container"
   echo ""
   exit 1 # Exit script after printing help
}

call_function() {
	#call the cript passed in parameters, save it in logs
	#end the program if the called script has returned an error
	
	function_command=$@
	log_script_call "$container_path" "$function_command"
	$function_command #execute the command
	if [ ! $? = 0 ] #test if command failed
	then
		echo "error in $(basename $2)"
		echo "canceling dna_store"
		exit 1
	fi
}


#----- default parameters -----#
n_synth=100
i_error=0
d_error=0
s_error=0

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

while [[ "$#" -gt 0 ]]; do
	case "$1" in
	    -n_synth ) if [ "$#" -lt 2 ] ; then
	    				echo "no value after param n_synth" ; exit 1 ;
	    			else
	    				n_synth="$2" ; shift 2 ;
	    			fi ;;
	    -i_error ) if [ "$#" -lt 2 ] ; then
	    				echo "no value after param i_error" ; exit 1 ;
	    			else
	    				i_error="$2" ; shift 2 ;
	    			fi ;;
	    -d_error ) if [ "$#" -lt 2 ] ; then
	    				echo "no value after param d_error" ; exit 1 ;
	    			else
	    				d_error="$2" ; shift 2 ;
	    			fi ;;
	    -s_error ) if [ "$#" -lt 2 ] ; then
	    				echo "no value after param s_error" ; exit 1 ;
	    			else
	    				s_error="$2" ; shift 2 ;
	    			fi ;;
	    -h | --help ) help_function ; exit 1;;
	    -* ) echo "unknown parameter $1" ; exit 1;;
	    * ) if test -z "${container_path:-}" ; then
				container_path="${1}" ; shift ;
			else
				echo "unknown parameter $1" ; exit 1 ;
			fi ;;
	esac
done

int_regex='^[0-9]+$'
if ! [[ $n_synth =~ $int_regex ]] ; then
   echo "error: n_synth ($n_synth) is not a number" ; exit 1
fi

float_regex='^[-+]?[0-9]+\.?[0-9]*$'
if ! [[ $i_error =~ $float_regex ]] ; then
   echo "error: i_error ($i_error) is not a float" ; exit 1
fi
float_regex='^[-+]?[0-9]+\.?[0-9]*$'
if ! [[ $d_error =~ $float_regex ]] ; then
   echo "error: d_error ($d_error) is not a float" ; exit 1
fi
float_regex='^[-+]?[0-9]+\.?[0-9]*$'
if ! [[ $s_error =~ $float_regex ]] ; then
   echo "error: s_error ($s_error) is not a float" ; exit 1
fi

if test -z "${container_path:-}"
then
	echo "container name missing"
	help_function
	exit 1
fi

if [ ! -d "$container_path" ] 
then
    echo "error : container $container_path not found" 
    exit 1
fi


#----------------------------------------------#
######### ====== store document ====== #########
#----------------------------------------------#
time=$(date +"%s")

is_editable=$(get_container_param $container_path "editable")

if ! $is_editable
then
	echo "the container is not editable"
	exit 1
fi


#----Fragments Design----#

#copy all the fragments from all the documents into one file
#cat "$container_path"/*/3_final_fragments.fasta > "$container_path"/container_fragments.fasta

# get the fragment length from the container options
frag_length=$(get_container_param "$container_path" "frag_length")

fragment_design_script="$project_dir"/synthesis_modules/ordered_fragment_design.py

for directory in "$container_path"/*/ ; do
	call_function python3 "$fragment_design_script" -i "$directory"/1_channel.fasta -o "$directory"/2_fragments.fasta -l $frag_length --cont_doc "$directory"
done

#----Synthesis----#

simulation=$(get_container_param $container_path "simulation")

if [ $simulation ]
then
	synthesis_script="$project_dir"/synthesis_modules/synthesis_simulation.py
	
	for directory in "$container_path"/*/ ; do
		call_function python3 "$synthesis_script" -i "$directory"/2_fragments.fasta -o "$directory"/3_synthesis.fasta -n $n_synth --i_error $i_error --d_error $d_error --s_error $s_error
	done
else
	echo "container is not in simulation mode"
	exit 0
	#TODO call a real synthesis
fi


#----Fragment assembly----#

if [ $simulation ]
then
	molecule_design_script="$project_dir"/synthesis_modules/assembly_simulation.py
	
	assembly_type=$(get_container_param $container_path "assembly_type")
	
	for directory in "$container_path"/*/ ; do
		document_index=$(basename $directory)
		n_frag=$(get_doc_param $container_path $document_index "fragment_number")
		n_mol=$(($n_frag * 200)) #TODO
				
		call_function python3 "$molecule_design_script" -i "$directory"/3_synthesis.fasta -o "$directory"/4_molecules.fasta -t $assembly_type -f $n_frag -n $n_mol
	done
	#concatenate all the molecules files into one to represent the physical container
	cat "$container_path"/*/4_molecules.fasta > "$container_path"/container_molecules.fasta
else
	echo "container is not in simulation mode"
	exit 0
fi


#----Update metadata----#

set_container_param $container_path "editable" false

echo "Documents of $(basename $container_path) successfully stored !"
end_time=$(date +"%s")
for directory in "$container_path"/*/ ; do
	echo "dna_store : $(($end_time - $time)) s" >> "$directory"/workflow_times.txt
done
exit 0
