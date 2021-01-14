#!/bin/bash

#--------------------------------------------------------------#
######### ===== Set the parameters in this part ====== #########
#--------------------------------------------------------------#

process_name="workflow_2_perfect_test" #name of the directory to save the generated files

#----- parameters for sequence generation -----#
random_seq=true #generate random sequences (true) or use an existing fasta file (false)
seq_path="workflow_2_fail/1_base_seq_file.fasta" #if random_seq is false, the sequences from this path are used (one .fasta file)
#else the sequences are generated with the following parameters
nbr_seq=1 #number of sequences
size_seq=1000 #size of the sequences
h_max=3 #maximum size for the homopolymeres

#----- parameters for overlap fragmentation -----#
frag_size=200 #size of the fragments
overlap_size=20 #size of the overlap

#----- parameters for synthesis -----#
spacer_path="spacer.fasta" #path to the spacer to se (.fasta file)
nbr_synth=500 #nomber of molecule to generate
n_frag=10 #number of sequence fragments in each molecule
i_error=0.00 #insertion error rate
d_error=0.00 #deletion error rate
s_error=0.00 #substitution error rate

#----- parameters for sequencing -----#
nbr_read=500 #number of read
perfect=2 # 0 = normal sequencing, 1 = no length repeat and noise, 2 = almost perfect reads without any randomness 

#-----------------------------------------------------#
######### ===== Part 0: initialisation ====== #########
#-----------------------------------------------------#
working_dir=$(pwd)
echo $HOME
if [ $HOME == "/Users/oboulle" ]
then
  project_dir="/Users/oboulle/Documents"
  conda_env="/Users/oboulle/anaconda2"
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
#start_time=$(($(date +%s%N)/1000000))

base_seq_path="$process_path/1_base_seq_file.fasta"

echo "Lancement du Workflow 2 : (taille $size_seq)"
echo "___Initialisation de la séquence de base___"

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
#end_time=$(($(date +%s%N)/1000000))
#echo "base sequences : $(($end_time - $start_time)) ms" >> $time_file
#------------------------------------------------------------#
######### ===== Part 2: overlap fragmentation ====== #########
#------------------------------------------------------------#
#start_time=$(($(date +%s%N)/1000000))
echo "___Fragmentation des séquences avec chevauchements___"

fragmentation_script="$project_dir/synthesis_simulation/overlap_fragmentation/overlap_fragmentation.py" #script for the overlap fragmentation
fragmentation_seq_path="$process_path/2_fragmented_seq_file.fasta"
python3 $fragmentation_script "$base_seq_path" "$fragmentation_seq_path" $frag_size $overlap_size
if [ ! $? = 0 ]
then
	exit 1
fi
#end_time=$(($(date +%s%N)/1000000))
#echo "overlap fragmentation : $(($end_time - $start_time)) ms" >> $time_file
#------------------------------------------------#
######### ===== Part 3: synthesis ====== #########
#------------------------------------------------#
#start_time=$(($(date +%s%N)/1000000))
echo "___Synthèse des séquences___"

synthesis_script="$project_dir/synthesis_simulation/synthesis_with_spacers/synthesis_with_spacers.py" #script for the synthesis with spacers
synthesis_path="$process_path/3_synthesis_file.fasta"
python3 $synthesis_script -i "$fragmentation_seq_path" -o "$synthesis_path" -s "$spacer_path" -n $nbr_synth --n_frag $n_frag --i_error $i_error --d_error $d_error --s_error $s_error
if [ ! $? = 0 ]
then
	exit 1
fi
#end_time=$(($(date +%s%N)/1000000))
#echo "synthesis : $(($end_time - $start_time)) ms" >> $time_file
#-----------------------------------------------------#
######### ===== Part 4: convert to fastq ====== #########
#-----------------------------------------------------#
fasta_to_fastq_script="$project_dir/synthesis_simulation/utils/dna_file_reader.py" #script to convert fasta to fastq
fastq_file="$process_path/4_sequencing.fastq"
python3 $fasta_to_fastq_script "$synthesis_path" "$fastq_file"
#-----------------------------------------------------#
######### ===== Part 5: reconstruction ====== #########
#-----------------------------------------------------#
#start_time=$(($(date +%s%N)/1000000))
echo "___Reconstruction___"
reconstruction_script="$project_dir/sequencing_simulation/spacer_sequencing/reconstruct.py" #script for the reconstruction
reconstruction_path="$process_path/5_result_sequence.fasta"
python3 $reconstruction_script "$fastq_file" "$reconstruction_path" "$spacer_path" $frag_size $overlap_size
if [ ! $? = 0 ]
then
	exit 1
fi
#end_time=$(($(date +%s%N)/1000000))
#echo "reconstruction : $(($end_time - $start_time)) ms" >> $time_file
#------------------------------------------------------#
######### ===== Part 6: Result Analysis ====== #########
#------------------------------------------------------#
echo "___Results___"
result_analysis_script="$project_dir/workflow_global/result_analysis/result_analysis_workflow_2.py"
python3 $result_analysis_script $base_seq_path $reconstruction_path

#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
