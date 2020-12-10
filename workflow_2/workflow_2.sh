#!/bin/bash

#--------------------------------------------------------------#
######### ===== Set the parameters in this part ====== #########
#--------------------------------------------------------------#

process_name="workflow_2_test" #name of the directory to save the generated files

#----- parameters for sequence generation -----#
random_seq=true #generate random sequences (true) or use an existing fasta file (false)
seq_path="test/1_base_seq_file.fasta" #if random_seq is false, the sequences from this path are used (one .fasta file)
#else the sequences are generated with the following parameters
nbr_seq=1 #number of sequences
size_seq=2000 #size of the sequences
h_max=3 #maximum size for the homopolymeres

#----- parameters for overlap fragmentation -----#
frag_size=220 #size of the fragments
overlap_size=20 #size of the overlap

#----- parameters for synthesis -----#
spacer_path="spacer.fasta" #path to the spacer to se (.fasta file)
nbr_synth=1000 #nomber of molecule to generate
n_frag=5 #number of sequence fragments in each molecule
i_error=0.1 #insertion error rate
d_error=0.00 #deletion error rate
s_error=0.00 #substitution error rate

#----- parameters for sequencing -----#
nbr_read=1000 #number of read
perfect=2 # 0 = normal sequencing, 1 = no length repeat and noise, 2 = almost perfect reads without any randomness 

#----- parameters for demultiplexing -----#
kmer_size=10 #size of the subsequences of the primers to search in the fastq sequences
point_threshold=10 #threshold of the minimum number of kmer found in the fastq_sequences to link it to a primer

#-----------------------------------------------------#
######### ===== Part 0: initialisation ====== #########
#-----------------------------------------------------#
working_dir=$(pwd)
if [ $HOME == "/home/oboulle" ]
then
  project_dir="/home/oboulle/Documents"
  conda_env="/home/oboulle/anaconda2"
else
  project_dir="/home/genouest/genscale/oboulle/documents"
	conda_env="/home/genouest/genscale/oboulle/anaconda2"
fi

process_path="$working_dir/$process_name"
rm -rf "$process_path"
mkdir -p "$process_path"
time_file="$process_path/times.txt"

#-----------------------------------------------------#
######### ===== Part 1: base sequences ====== #########
#-----------------------------------------------------#
start_time=$(($(date +%s%N)/1000000))

base_seq_path="$process_path/1_base_seq_file.fasta"

echo "Lancement du Workflow 2 :"
echo "___Initialisation des séquences de base___"

seq_gen_script="$project_dir/synthesis_simulation/sequence_generator/sequence_generator.py" #script for the sequence generation
if [ $random_seq = true ]
then
	python3 $seq_gen_script "$base_seq_path" $nbr_seq $size_seq $h_max
	if [ ! $? = 0 ]
	then
		exit 1
	fi
else
	if [ -f "$seq_path" ]
	then
		cp "$seq_path" "$base_seq_path"
	else
		echo "sequence file $seq_path not found !"
		exit 1
	fi
fi
end_time=$(($(date +%s%N)/1000000))
echo "base sequences : $(($end_time - $start_time)) ms" >> $time_file
#------------------------------------------------------------#
######### ===== Part 2: overlap fragmentation ====== #########
#------------------------------------------------------------#
start_time=$(($(date +%s%N)/1000000))
echo "___Fragmentation des séquences avec chevauchements___"

fragmentation_script="$project_dir/synthesis_simulation/overlap_fragmentation/overlap_fragmentation.py" #script for the overlap fragmentation
fragmentation_seq_path="$process_path/2_fragmented_seq_file.fasta"
python3 $fragmentation_script "$base_seq_path" "$fragmentation_seq_path" $frag_size $overlap_size
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(($(date +%s%N)/1000000))
echo "overlap fragmentation : $(($end_time - $start_time)) ms" >> $time_file
#------------------------------------------------#
######### ===== Part 3: synthesis ====== #########
#------------------------------------------------#
start_time=$(($(date +%s%N)/1000000))
echo "___Synthèse des séquences___"

synthesis_script="$project_dir/synthesis_simulation/synthesis_with_spacers/synthesis_with_spacers.py" #script for the synthesis with spacers
synthesis_path="$process_path/3_synthesis_file.fasta"
python3 $synthesis_script -i "$fragmentation_seq_path" -o "$synthesis_path" -s "$spacer_path" -n $nbr_synth --n_frag $n_frag --i_error $i_error --d_error $d_error --s_error $s_error
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(($(date +%s%N)/1000000))
echo "synthesis : $(($end_time - $start_time)) ms" >> $time_file
#-------------------------------------------------#
######### ===== Part 4: sequencing ====== #########
#-------------------------------------------------#
start_time=$(($(date +%s%N)/1000000))
echo "___Séquençage___"

deep_simu_home="$project_dir/sequencing_simulation/deep_simulator" #home of DeepSimulator
deep_simulator_script="$project_dir/sequencing_simulation/deep_simulator/deep_simulator.sh" #script for the sequencing
sequencing_path="$process_path/4_sequencing"
conda config --add envs_dirs "$conda_env/envs"
$deep_simulator_script -i "$synthesis_path" -o "$sequencing_path" -H "$deep_simu_home" -C "$conda_env" -n $nbr_read -P $perfect
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(($(date +%s%N)/1000000))
echo "sequencing : $(($end_time - $start_time)) ms" >> $time_file
#--------------------------------------------------#
######### ===== Part 5: basecalling ====== #########
#--------------------------------------------------#
start_time=$(($(date +%s%N)/1000000))
echo "___Base Calling___"

basecaller_path="$project_dir/sequencing_simulation/ont-guppy-cpu_4.2.2_linux64/ont-guppy-cpu/bin/guppy_basecaller" #path of the basecaller to use
basecalling_path="$process_path/5_basecalling"
mkdir "$basecalling_path"
$basecaller_path -r --input_path "$sequencing_path" --save_path "$basecalling_path" -c dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 8 --num_callers 1
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(($(date +%s%N)/1000000))
echo "basecalling : $(($end_time - $start_time)) ms" >> $time_file
#-----------------------------------------------------#
######### ===== Part 6: reconstruction ====== #########
#-----------------------------------------------------#
start_time=$(($(date +%s%N)/1000000))
echo "___Reconstruction___"
reconstruction_script="$project_dir/sequencing_simulation/spacer_sequencing/reconstruct.py" #script for the reconstruction
reconstruction_path="$process_path/6_reconstruction"
mkdir "$reconstruction_path"
fastq_file=(*.fastq)
python3 $reconstruction_script "$basecalling_path/"$fastq_file "$reconstruction_path/result_sequence.fasta" "$spacer_path" $frag_size $overlap_size
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(($(date +%s%N)/1000000))
echo "reconstruction : $(($end_time - $start_time)) ms" >> $time_file
#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
