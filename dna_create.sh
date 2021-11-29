#!/bin/bash

#DNA CREATE - initialise an empty container

help_function() {
   echo ""
   echo "Usage: dna_create [-sim] [-nocd] [-fl int] [-sp string] Cname"
   echo -e "\t-sim : turn on simulator mode"
   echo -e "\t-nocd : turn off channel encoding"
   echo -e "\t-fl : specify fragment length [default = 100]"
   echo -e "\t-sp : specify spacer [default = AAAAAAAACCCCCCCC]"
   echo -e "\tCname : name of the container"
   echo ""
   exit 1 # Exit script after printing help
}

#----- default parameters -----#
simulation=false
channel_coding=true
frag_length=100
spacer="AAAAAAAACCCCCCCC"

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

while true; do
  case "$1" in
    -sim ) simulation=true ; shift ;;
    -nocd ) channel_coding=false ; shift ;;
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

echo "sim=$simulation channel_coding=$channel_coding frag_length=$frag_length spacer=$spacer container=$container_name"

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

options_file="$container_path"/.options

cat > "$options_file" << eof
simulation $simulation
channel_coding $channel_coding
frag_length $frag_length
spacer $spacer
container_name $container_name
eof

cdi_file="$container_path"/.cdi

cat > "$cdi_file" << eof
0
eof

echo "Container $container_name created successfully !"

exit 0
