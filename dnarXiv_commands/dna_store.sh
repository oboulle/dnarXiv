#!/bin/bash

#----- default parameters -----#
n_synth=100
i_error=0
d_error=0
s_error=0

#-----------------------------------------------#
######### ====== read parameters ====== #########
#-----------------------------------------------#

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

while true; do
  case "$1" in
    -n_synth ) n_synth="$2" ; shift 2 ;;
    -i_error ) i_error="$2" ; shift 2 ;;
    -d_error ) d_error="$2" ; shift 2 ;;
    -s_error ) s_error="$2" ; shift 2 ;;
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

#--------------------------------------------#
######### ====== store document ====== #########
#--------------------------------------------#

cdi_file="$container_path/.cdi"
container_index=$(head -n 1 "$cdi_file")

if (( $container_index < 0 ))
then
	echo "the container is not editable"
	exit 1
fi

if [ $HOME == "/Users/oboulle" ]
then
  project_dir="/Users/oboulle/Documents"
else
  project_dir="/home/oboulle/Documents"
fi

cancel_dna_store() {
	#TODO
	#remove primers
	#remove META dir
	echo "cancel dna_store"
}

# get parameters from the container options
while read var value; do
    export "$var"="$value"
done < "$container_path/.options"

#----Primers Generation----#

primers_generation_script="$project_dir/synthesis_simulation/primer_generation/primer_generation.py"

for directory in $container_path/*/ ; do
	python3 $primers_generation_script "$directory/fragments.fasta" "$directory/primers.fasta" "$spacer" "$directory/.meta"
	if [ ! $? = 0 ]
	then
		echo "error in primers generation"
		#cancel_dna_store#TODO
		#exit 1
	fi
done

#----Global metadata generation----#

meta_doc_path="$container_path/META"

if [ -d "$meta_doc_path" ] 
then
    rm -rf $meta_doc_path
fi 

#create a META dir in the container
mkdir -p "$meta_doc_path"

metadata_generation_script="$project_dir/synthesis_simulation/source_encoding/metadata_generation.py" 

meta_concatenation_path="$meta_doc_path/concatenation.txt"

#generate the concatenation of the container .meta files 
python3 $metadata_generation_script "$container_path" "$meta_concatenation_path" #TODO
if [ ! $? = 0 ]
then
	echo "error in global metadata generation"
	#cancel_dna_store#TODO
	#exit 1
fi

#add the metadata concatenation to the container as a file (source encoding, channel encoding, etc)
$project_dir/workflow_global/dnarXiv_commands/dna_add.sh -meta "$meta_concatenation_path" "$container_path"
if [ ! $? = 0 ]
then
	echo "error in global metadata addition"
	#cancel_dna_store#TODO
	#exit 1
fi

#generate the primers for this new file
python3 $primers_generation_script "$meta_doc_path/source.fasta" "$meta_doc_path/primers.fasta" "$spacer" "$meta_doc_path/.meta"
if [ ! $? = 0 ]
then
	echo "error in primers generation"
	#cancel_dna_store#TODO
	#exit 1
fi

#----Synthesis----#

if [ "$simulation" = true ]
then
	synthesis_script="$project_dir/synthesis_simulation/synthesis/synthesis.py"
	
	for directory in $container_path/*/ ; do
		python3 $synthesis_script -i "$directory/fragments.fasta" -o "$directory/synthesis.fasta" -n $n_synth --i_error $i_error --d_error $d_error --s_error $s_error
		if [ ! $? = 0 ]
		then
			echo "error in synthesis simulation"
			#cancel_dna_store#TODO
			#exit 1
		fi
	done
else
	echo "container is not in simulation mode"
	#TODO call a real synthesis
fi

#----Molecule design----#

if [ "$simulation" = true ]
then
	molecule_design_script="$project_dir/synthesis_simulation/synthesis/molecule_design.py"
	
	for directory in $container_path/*/ ; do
		python3 $molecule_design_script -i "$directory/synthesis.fasta" -o "$directory/molecules.fasta" -s "$spacer" -p "$directory/primers.fasta" -n $n_synth
		if [ ! $? = 0 ]
		then
			echo "error in molecule design"
			#cancel_dna_store#TODO
			#exit 1
		fi
	done
else
	echo ""
fi

#----Update .cdi----#

cat > $cdi_file << eof
-1
eof

echo "Documents of $container_path successfully stored !"
exit 0
