#!/bin/bash

#--------------------------------------------------------------#
######### ===== Set the parameters in this part ====== #########
#--------------------------------------------------------------#
project_dir="/home/genouest/genscale/oboulle/documents"
process_name="shell_script_test22" #name of the directory to save the generated files
workflow_path="$project_dir/workflow_global/workflow_1" #path of the workflow directory

#----- parameters for sequence generation -----#
random_seq=true #generate random sequences (true) or use an existing fasta file (false)
seq_path="test/1_base_seq_file.fasta" #if random_seq is false, the sequences from this path are used (one .fasta file)
#else the sequences are generated
seq_gen_script="$project_dir/synthesis_simulation/sequence_generator/sequence_generator.py" #script for the sequence generation
nbr_seq="3" #number of sequences
size_seq="200" #size of the sequences
h_max="3" #maximum size for the homopolymeres

#----- parameters for primer addition -----#
primer_script="$project_dir/synthesis_simulation/sequence_primer/sequence_primer.py" #script for the primer addition
#primer_path="$project_dir/synthesis_simulation/sequence_primer/test_primer.fasta" #path of the primer file (.fasta)

#----- parameters for synthesis -----#
synthesis_script="$project_dir/synthesis_simulation/synthesis_simulator.py" #script for the synthesis
nbr_synth=1000 #nomber of synthesis for each sequence
i_error=0.00 #insertion error rate
d_error=0.00 #deletion error rate
s_error=0.00 #substitution error rate

#----- parameters for sequencing -----#
deep_simu_home="$project_dir/sequencing_simulation/deep_simulator" #home of DeepSimulator
deep_simulator_script="$project_dir/sequencing_simulation/deep_simulator/deep_simulator.sh" #script for the sequencing
conda_env="/home/genouest/genscale/oboulle/anaconda2"
nbr_read=100 #number of read
perfect=2 # 0 = normal sequencing, 1 = no length repeat and noise, 2 = almost perfect reads without any randomness 

#----- parameters for basecalling -----#
basecaller_path="$project_dir/sequencing_simulation/ont-guppy-cpu_4.2.2_linux64/ont-guppy-cpu/bin/guppy_basecaller" #path of the basecaller to use

#----- parameters for demultiplexing -----#
demultiplexing_script="$project_dir/sequencing_simulation/demultiplexing/demultiplexing.py" #script for the demultiplexing
kmer_size=10 #size of the subsequences of the primers to search in the fastq sequences
point_threshold=10 #threshold of the minimum number of kmer found in the fastq_sequences to link it to a primer

#----- parameters for consensus -----#
consensus_script="$project_dir/sequencing_post_processing/ccsa.py"

#-----------------------------------------------------#
######### ===== Part 1: base sequences ====== #########
#-----------------------------------------------------#
working_dir=$(pwd)

process_path="$working_dir/$process_name"
echo $process_path
rm -rf $process_path
mkdir -p $process_path

base_seq_path="$process_path/1_base_seq_file.fasta"

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

primer_seq_path="$process_path/2_primer_seq_file.fasta"
python3 $primer_script $base_seq_path $primer_seq_path $process_path
if [ ! $? = 0 ]
then
	exit 1
fi

#------------------------------------------------#
######### ===== Part 3: synthesis ====== #########
#------------------------------------------------#

echo "___Synthèse des séquences___"

synthesis_path="$process_path/3_synthesis_file.fasta"
python3 $synthesis_script -i $primer_seq_path -o $synthesis_path -n $nbr_synth --i_error $i_error --d_error $d_error --s_error $s_error
if [ ! $? = 0 ]
then
	exit 1
fi

#-------------------------------------------------#
######### ===== Part 4: sequencing ====== #########
#-------------------------------------------------#

echo "___Séquençage___"

sequencing_path="$process_path/4_sequencing"
conda config --add envs_dirs "$conda_env/envs"
$deep_simulator_script -i $synthesis_path -o $sequencing_path -H $deep_simu_home -C $conda_env -n $nbr_read -P $perfect
if [ ! $? = 0 ]
then
	exit 1
fi
	
#--------------------------------------------------#
######### ===== Part 5: basecalling ====== #########
#--------------------------------------------------#

echo "___Base Calling___"

basecalling_path="$process_path/5_basecalling"
mkdir $basecalling_path
$basecaller_path -r --input_path $sequencing_path --save_path $basecalling_path -c dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 8 --num_callers 1
if [ ! $? = 0 ]
then
	exit 1
fi

#-----------------------------------------------------#
######### ===== Part 6: demultiplexing ====== #########
#-----------------------------------------------------#

echo "___Demultiplexing___"

demultiplexing_path="$process_path/6_demultiplexing"
fastq_file=(*.fastq)
python3 $demultiplexing_script -i $basecalling_path/$fastq_file -o $demultiplexing_path -p $process_path/primers --kmer_size $kmer_size --point-threshold $point_threshold
if [ ! $? = 0 ]
then
	exit 1
fi

#------------------------------------------------#
######### ===== Part 7: Consensus ====== #########
#------------------------------------------------#
echo "___Consensus___"

consensus_path="$process_path/7_consensus"
mkdir $consensus_path
for sequence_file in $demultiplexing_path/*.fastq
do
  sequence_name=`eval basename -s .fastq $sequence_file`
	if [ $sequence_name != unlinked ]
	then
	  python3 $consensus_script -read $sequence_file -primer $process_path/primers/$sequence_name.fasta -length $size_seq -out $consensus_path/$sequence_name.fasta -graph $consensus_path/$sequence_name.gexf
	fi
done

#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
