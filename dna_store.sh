#!/bin/bash
set -u #exit and display error message if a variable is empty 

#DNA STORE - synthetise the documents of the container into molecules

#set -x #display each command used

help_function() {
   echo ""
   echo "Usage: dna_store [-n_synth int] [-i_error int] [-d_error int] [-s_error int] Cname"
   echo -e "\t-n_synth : specify number of synthesis per fragments [default = 100]"
   echo -e "\t-i_error : specify insertion error rate in synthesis simulation [default = 0]"
   echo -e "\t-d_error : specify deletion error rate in synthesis simulation [default = 0]"
   echo -e "\t-s_error : specify substitution error rate in synthesis simulation [default = 0]"
   echo -e "\tCname : path to the container"
   echo ""
   exit 1 # Exit script after printing help
}

check_error_function () { #end the program if the previously called script has returned an error
	if [ ! $? = 0 ]
	then
		echo "error in $1"
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

while true; do
  case "$1" in
    -n_synth ) n_synth="$2" ; shift 2 ;;
    -i_error ) i_error="$2" ; shift 2 ;;
    -d_error ) d_error="$2" ; shift 2 ;;
    -s_error ) s_error="$2" ; shift 2 ;;
    -h | --help ) help_function ; exit 1;;
    -* ) echo "unknown parameter $1" ; exit 1;;
    * ) container_path="${1}" ; break ;;
  esac
done

int_regex='^[0-9]+$'
if ! [[ $n_synth =~ $int_regex ]] ; then
   echo "error: n_synth ($n_synth) is not a number" ; exit 1
fi

if test -z "$container_path"
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

project_dir="$(dirname $0)/.." #parent of the directory containing this script

source "$project_dir"/workflow_commands/metadata_manager.sh #load the xml manager script
meta_file="$container_path"/metadata.xml

is_editable=$(get_container_param $meta_file "editable")

if ! $is_editable
then
	echo "the container is not editable"
	exit 1
fi


#----Primers Generation----#

#copy all the fragments from all the documents into one file
cat "$container_path"/*/3_final_fragments.fasta > "$container_path"/container_fragments.fasta

primers_generation_script="$project_dir"/synthesis_simulation/primer_generation.py

python3 "$primers_generation_script" "$container_path"
check_error_function "primer generation"

simulation=$(get_container_param $meta_file "simulation")


#----Synthesis----#

if [ $simulation ]
then
	synthesis_script="$project_dir"/synthesis_simulation/synthesis.py
	
	for directory in "$container_path"/*/ ; do
		python3 "$synthesis_script" -i "$directory"/3_final_fragments.fasta -o "$directory"/4_synthesis.fasta -n $n_synth --i_error $i_error --d_error $d_error --s_error $s_error
		check_error_function "synthesis simulation"
	done
else
	echo "container is not in simulation mode"
	exit 0
	#TODO call a real synthesis
fi


#----Molecule design----#

if [ $simulation ]
then
	molecule_design_script="$project_dir"/synthesis_simulation/molecule_design.py
	for directory in "$container_path"/*/ ; do
		document_index=$(basename $directory)
		n_frag=$(get_doc_param $meta_file $document_index "fragment_number")
		start_primer=$(get_doc_param $meta_file $document_index "start_primer")
		stop_primer=$(get_doc_param $meta_file $document_index "stop_primer")
		n_mol=$(($n_frag * 1000)) #TODO
				
		python3 "$molecule_design_script" -i "$directory"/4_synthesis.fasta -o "$directory"/5_molecules.fasta --start $start_primer --stop $stop_primer -n $n_mol
		check_error_function "error in molecule design"
	done
	#concatenate all the molecules files into one to represent the physical container
	cat "$container_path"/*/5_molecules.fasta > "$container_path"/container_molecules.fasta
else
	echo "container is not in simulation mode"
	exit 0
fi


#----Update metadata----#

set_container_param $meta_file "editable" false

echo "Documents of $container_path successfully stored !"
end_time=$(date +"%s")
for directory in "$container_path"/*/ ; do
	echo "dna_store : $(($end_time - $time)) s" >> "$directory"/workflow_times.txt
done
exit 0
