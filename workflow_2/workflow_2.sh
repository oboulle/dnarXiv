#!/bin/bash

#--------------------------------------------------------------#
######### ===== Set the parameters in this part ====== #########
#--------------------------------------------------------------#

process_name="workflow_2_test" #name of the directory to save the generated files

#----- parameters for sequence generation -----#
random_seq=true #generate random sequences (true) or use an existing fasta file (false)
seq_path="workflow_2_fail/1_base_seq_file.fasta" #if random_seq is false, the sequences from this path are used (one .fasta file)
#else the sequences are generated with the following parameters
nbr_seq=1 #number of sequences
size_seq=10100 #size of the sequences
h_max=3 #maximum size for the homopolymeres

#----- parameters for overlap fragmentation -----#
frag_size=200 #size of the fragments
overlap_size=20 #size of the overlap

#----- parameters for synthesis -----#
spacer_path="spacer.fasta" #path to the spacer to se (.fasta file)
nbr_synth=300 #nomber of molecule to generate
n_frag=10 #number of sequence fragments in each molecule
i_error=0.00 #insertion error rate
d_error=0.00 #deletion error rate
s_error=0.00 #substitution error rate

#----- parameters for sequencing -----#
nbr_read=300 #number of read
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
fi
if [ $HOME == "/home/oboulle" ]
then
  project_dir="/home/oboulle/Documents"
  conda_env="/home/oboulle/anaconda2"
fi

process_path="$working_dir/$process_name"_$(date +"%H:%M:%S")
rm -rf "$process_path"
mkdir -p "$process_path"
time_file="$process_path/times.txt"

summary="$process_path/summary.txt"

cat > $summary << eof
$process_name

--parameters--
size_seq : $size_seq
frag_size : $frag_size
overlap_size : $overlap_size
nbr_synth : $nbr_synth
n_frag : $n_frag
i_error : $i_error
d_error : $d_error
s_error : $s_error
nbr_read : $nbr_read
perfect_sequencing : $perfect

--results--
eof

#-----------------------------------------------------#
######### ===== Part 1: base sequences ====== #########
#-----------------------------------------------------#
start_time=$(date +"%s")

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
end_time=$(date +"%s")
echo "base sequences : $(($end_time - $start_time)) s" >> $time_file
#------------------------------------------------------------#
######### ===== Part 2: overlap fragmentation ====== #########
#------------------------------------------------------------#
start_time=$(date +"%s")
echo "___Fragmentation des séquences avec chevauchements___"

fragmentation_script="$project_dir/synthesis_simulation/overlap_fragmentation/overlap_fragmentation.py" #script for the overlap fragmentation
fragmentation_seq_path="$process_path/2_fragmented_seq_file.fasta"
python3 $fragmentation_script "$base_seq_path" "$fragmentation_seq_path" $frag_size $overlap_size
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(date +"%s")
echo "overlap fragmentation : $(($end_time - $start_time)) ms" >> $time_file
#------------------------------------------------#
######### ===== Part 3: synthesis ====== #########
#------------------------------------------------#
start_time=$(date +"%s")
echo "___Synthèse des séquences___"

synthesis_script="$project_dir/synthesis_simulation/synthesis_with_spacers/synthesis_with_spacers.py" #script for the synthesis with spacers
synthesis_path="$process_path/3_synthesis_file.fasta"
python3 $synthesis_script -i "$fragmentation_seq_path" -o "$synthesis_path" -s "$spacer_path" -n $nbr_synth --n_frag $n_frag --i_error $i_error --d_error $d_error --s_error $s_error
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(date +"%s")
echo "synthesis : $(($end_time - $start_time)) s" >> $time_file
#-------------------------------------------------#
######### ===== Part 4: sequencing ====== #########
#-------------------------------------------------#
start_time=$(date +"%s")
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
end_time=$(date +"%s")
echo "sequencing : $(($end_time - $start_time)) s" >> $time_file
#--------------------------------------------------#
######### ===== Part 5: basecalling ====== #########
#--------------------------------------------------#
start_time=$(date +"%s")
echo "___Base Calling___"

basecaller_path="$project_dir/sequencing_simulation/ont-guppy-cpu/bin/guppy_basecaller" #path of the basecaller to use
basecalling_path="$process_path/5_basecalling"
mkdir "$basecalling_path"
$basecaller_path -r --input_path "$sequencing_path" --save_path "$basecalling_path" -c dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 8 --num_callers 1
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(date +"%s")
echo "basecalling : $(($end_time - $start_time)) s" >> $time_file
#-----------------------------------------------------#
######### ===== Part 6: reconstruction ====== #########
#-----------------------------------------------------#
start_time=$(date +"%s")
echo "___Reconstruction___"
reconstruction_script="$project_dir/sequencing_simulation/spacer_sequencing/reconstruct.py" #script for the reconstruction
reconstruction_path="$process_path/6_result_sequence.fasta"
fastq_file=(*.fastq)
python3 $reconstruction_script "$basecalling_path/"$fastq_file "$reconstruction_path" "$spacer_path" $frag_size $overlap_size
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(date +"%s")
echo "reconstruction : $(($end_time - $start_time)) s" >> $time_file
#------------------------------------------------------#
######### ===== Part 7: Result Analysis ====== #########
#------------------------------------------------------#
echo "___Results___"
result_analysis_script="$project_dir/workflow_global/result_analysis/result_analysis_workflow_2.py"
python3 $result_analysis_script $base_seq_path $reconstruction_path | while read line ; do
    echo $line >> $summary
done

#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
