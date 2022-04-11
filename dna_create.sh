#!/bin/bash
set -u #exit and display error message if a variable is empty 

#DNA CREATE - initialise an empty container

project_dir="$(dirname ${BASH_SOURCE})/.." #parent of the directory containing this script
source "$project_dir"/workflow_commands/metadata_manager.sh #load the xml manager script

help_function() {
   echo ""
   echo "Usage: dna_create [-sim] [-fl int] [-sp string] Cname"
   echo -e "\t-sim : turn on simulator mode"
   echo -e "\t-fl : specify fragment length [default = 100]"
   echo -e "\t-sp : specify spacer [default = AAAAAAAACCCCCCCC]"
   echo -e "\tCname : name of the container"
   echo ""
   exit 1 # Exit script after printing help
}

#----- default parameters -----#
simulation=true
frag_length=100
spacer="AAAAAAAACCCCCCCC"

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

while true; do
  case "$1" in
    -sim ) simulation=true ; shift ;;
    -fl ) frag_length=$2 ; shift 2 ;;
    -sp ) spacer=$2 ; shift 2 ;;
    -h | --help ) help_function ; exit 1;;
    -* ) echo "unknown parameter $1" ; exit 1;;
    * ) container_name="${1}" ; break ;;
  esac
done

int_regex='^[0-9]+$'
if ! [[ $frag_length =~ $int_regex ]]
then
   echo "error: fragment length ($frag_length) is not a number" ; exit 1
fi

dna_regex='^[acgtACGT]+$'
if ! [[ $spacer =~ $dna_regex ]]
then
   echo "error: spacer ($spacer) is not a dna sequence" ; exit 1
fi

if test -z "$container_name"
then
	echo "container name missing"
	help_function
	exit 1
fi

echo "sim=$simulation frag_length=$frag_length spacer=$spacer container=$container_name"

#----------------------------------------------------------#
######### ====== create directory and files ====== #########
#----------------------------------------------------------#
container_path="$(pwd)/$container_name"

if [ -d "$container_path" ] 
then
    echo "error : a container named $container_name already exists" 
    exit 1
fi 

mkdir -p "$container_path"

init_metadata_file "$container_path" "$container_name" $simulation $frag_length "ordered" false "$spacer" "nanopore"

echo "Container $container_name created successfully !"

exit 0
