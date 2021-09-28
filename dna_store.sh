#!/bin/bash

#DNA STORE - synthetise the documents of the container into molecules

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

cdi_file="$container_path"/.cdi
container_index=$(head -n 1 "$cdi_file")

if (( $container_index < 0 ))
then
	echo "the container is not editable"
	exit 1
fi

if [ $HOME == "/Users/oboulle" ]
then
	project_dir="/Users/oboulle/Documents"
elif [ $HOME == "/udd/oboulle" ]
then
	project_dir="/udd/oboulle/Documents"
else
	project_dir="/home/oboulle/Documents"
fi

# get parameters from the container options
while read var value; do
    export "$var"="$value"
done < "$container_path"/.options

#----Primers Generation----#

#copy all the fragments from all the documents into one file
cat "$container_path"/*/3_final_fragments.fasta > "$container_path"/container_fragments.fasta

primers_generation_script="$project_dir"/synthesis_simulation/primer_generation.py

python3 "$primers_generation_script" "$container_path"
check_error_function "primer generation"

#----Synthesis----#

if [ $simulation = true ]
then
	synthesis_script="$project_dir"/synthesis_simulation/synthesis.py
	
	for directory in "$container_path"/*/ ; do
		python3 "$synthesis_script" -i "$directory"/3_final_fragments.fasta -o "$directory"/4_synthesis.fasta -n $n_synth --i_error $i_error --d_error $d_error --s_error $s_error
		check_error_function "synthesis simulation"
	done
else
	echo "container is not in simulation mode"
	#TODO call a real synthesis
fi

#----Molecule design----#

if [ $simulation = true ]
then
	molecule_design_script="$project_dir"/synthesis_simulation/molecule_design.py
	for directory in "$container_path"/*/ ; do
		# get parameters from the .meta file of the document to store
		while read var value; do
		    export "$var"="$value"
		done < "$directory"/.meta
		n_mol=$(($n_frag * 50)) #TODO
		python3 "$molecule_design_script" -i "$directory"/4_synthesis.fasta -o "$directory"/5_molecules.fasta -s $spacer -p "$directory"/primers.fasta -n $n_mol
		check_error_function "error in molecule design"
	done
	#concatenate all the molecules files into one to represent the physical container
	cat "$container_path"/*/5_molecules.fasta > "$container_path"/container_molecules.fasta
else
	echo ""
fi

#----Update .cdi----#

cat > "$cdi_file" << eof
-1
eof

echo "Documents of $container_path successfully stored !"
end_time=$(date +"%s")
for directory in "$container_path"/*/ ; do
	echo "dna_store : $(($end_time - $time)) s" >> "$directory"/workflow_times.txt
done
exit 0
