#!/bin/bash

#--------------------------------------------------------------#
######### ===== Set the parameters in this part ====== #########
#--------------------------------------------------------------#

process_name="exemple" #name of the directory to save the generated files

#----- parameters for sequence generation -----#
random_seq=true #generate random sequences (true) or use an existing fasta file (false)
seq_path="test/1_base_seq_file.fasta" #if random_seq is false, the sequences from this path are used (one .fasta file)
#else the sequences are generated
seq_gen_script="/udd/oboulle/Documents/synthesis_simulation/sequence_generator/sequence_generator.py" #script for the sequence generation
nbr_seq="3" #number of sequences
size_seq="100" #size of the sequences
h_max="3" #maximum size for the homopolymeres

#----- parameters for primer addition -----#
primer_script="/udd/oboulle/Documents/synthesis_simulation/sequence_primer/sequence_primer.py" #script for the primer addition
primer_path="/udd/oboulle/Documents/synthesis_simulation/sequence_primer/test_primer.fasta" #path of the primer file (.fasta)

#----- parameters for synthesis -----#
synthesis_script="/udd/oboulle/Documents/synthesis_simulation/synthesis_simulator.py" #script for the synthesis
nbr_synth=1000 #nomber of synthesis for each sequence
i_error=0.01 #insertion error rate
d_error=0.01 #deletion error rate
s_error=0.01 #substitution error rate

#----- parameters for sequencing -----#
deep_simu_home="/udd/oboulle/Documents/sequencing_simulation/deep_simulator" #home of DeepSimulator
deep_simulator_script="/udd/oboulle/Documents/sequencing_simulation/deep_simulator/deep_simulator.sh" #script for the sequencing
nbr_read=100 #number of read

#----- parameters for basecalling -----#
basecaller_path="/udd/oboulle/Documents/sequencing_simulation/basecalling/guppy-cpu_4.2.2/bin/guppy_basecaller" #path of the basecaller to use

#-----------------------------------------------------#
######### ===== Part 1: base sequences ====== #########
#-----------------------------------------------------#

rm -rf $process_name
mkdir -p $process_name

base_seq_path="$process_name/1_base_seq_file.fasta"

echo "___Initialisation des séquences de base___"

if [ $random_seq = true ]
then
	python3 $seq_gen_script $base_seq_path $nbr_seq $size_seq $h_max
	if [ ! $? = 0 ]
	then
		exit 1
	fi
else
	if [ -f "$seq_path" ]
	then
		cp $seq_path $base_seq_path
	else
		echo "sequence file $seq_path not found !"
		exit 1
	fi
fi

#-------------------------------------------------------#
######### ===== Part 2: primers addition ====== #########
#-------------------------------------------------------#

echo "___Ajout de primers en début et fin de séquences___"

primer_seq_path="$process_name/2_primer_seq_file.fasta"
python3 $primer_script $base_seq_path $primer_seq_path $primer_path
if [ ! $? = 0 ]
then
	exit 1
fi

#------------------------------------------------#
######### ===== Part 3: synthesis ====== #########
#------------------------------------------------#

echo "___Synthèse des séquences___"

synthesis_path="$process_name/3_synthesis_file.fasta"
python3 $synthesis_script -i $primer_seq_path -o $synthesis_path -n $nbr_synth --i_error $i_error --d_error $d_error --s_error $s_error
if [ ! $? = 0 ]
then
	exit 1
fi

#-------------------------------------------------#
######### ===== Part 4: sequencing ====== #########
#-------------------------------------------------#

echo "___Séquençage___"

sequencing_path="$process_name/4_sequencing"
$deep_simulator_script -i $synthesis_path -o $sequencing_path -H $deep_simu_home -n $nbr_read
if [ ! $? = 0 ]
then
	exit 1
fi
	
#--------------------------------------------------#
######### ===== Part 5: basecalling ====== #########
#--------------------------------------------------#

echo "___Base Calling___"

basecalling_path="$process_name/5_basecalling"
mkdir $basecalling_path
$basecaller_path -r --input_path $sequencing_path --save_path $basecalling_path -c dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 8 --num_callers 1
if [ ! $? = 0 ]
then
	exit 1
fi
	
#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
