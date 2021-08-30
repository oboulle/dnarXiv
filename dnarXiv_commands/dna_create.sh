#!/bin/bash

help_function() {
   echo ""
   echo "Usage: dna_create [-sim] [-fl int] [-sp string] Cname"
   echo -e "\t-sim : turn on simulator mode [default = false]"
   echo -e "\t-fl : specify fragment length [default = 200]"
   echo -e "\t-sp : specify spacer [default = AAAAAAAACCCCCCCC]"
   echo -e "\tCname : name of the container"
   echo ""
   exit 1 # Exit script after printing help
}

#----- default parameters -----#
simulation=false
frag_length=200
spacer="AAAAAAAACCCCCCCC"

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

while true; do
  case "$1" in
    -sim ) simulation=true ; shift ;;
    -fl ) frag_length="$2" ; shift 2 ;;
    -sp ) spacer="$2" ; shift 2 ;;
    -h | --help ) help_function ; shift ;;
    -* ) echo "unknown parameter $1" ; exit 1;;
    * ) container_name="${1}" ; break ;;
  esac
done

int_regex='^[0-9]+$'
if ! [[ $frag_length =~ $int_regex ]] ; then
   echo "error: fragment length ($frag_length) is not a number" ; exit 1
fi

dna_regex='^[acgtACGT]+$'
if ! [[ $spacer =~ $dna_regex ]] ; then
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
echo workflow $(date +"%Hh%Mm%S") >> workflow_times.txt
time=$(date +"%s")

container_path="$(pwd)/$container_name"

if [ -d "$container_path" ] 
then
    echo "error : a container named $container_name already exists" 
    exit 1
fi 

mkdir -p "$container_path"

options_file="$container_path/.options"

cat > $options_file << eof
simulation $simulation
frag_length $frag_length
spacer $spacer
container_name $container_name
eof

cdi_file="$container_path/.cdi"

cat > $cdi_file << eof
0
eof

echo "Container $container_name created successfully !"
end_time=$(date +"%s")
echo "dna_create : $(($end_time - $time)) s" >> workflow_times.txt
exit 0
