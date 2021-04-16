#!/bin/bash

#--------------------------------------------------------------#
######### ===== Set the parameters in this part ====== #########
#--------------------------------------------------------------#

process_name="workflow_4_test/workflow_4" #name of the directory to save the generated files

#----- parameters for sequence generation -----#
random_seq=true #generate random sequences (true) or use an existing fasta file (false)
seq_path="workflow_4_fail/1_base_seq_file.fasta" #if random_seq is false, the sequences from this path are used (one .fasta file)
#else the sequences are generated with the following parameters
nbr_seq=1 #number of sequences
size_seq=1000 #size of the sequence
h_max=3 #maximum size for the homopolymeres

#----- parameters for fragmentation -----#
frag_size=200 #size of the fragments
tag_size=13 #size of the tag (4+3*n for 16^(n+1) maximum number)
spacer_path="spacer.fasta" #path to the spacer to use (.fasta file)

#----- parameters for synthesis -----#
primers_path="primers.fasta" #path to the primers to use (.fasta file)
nbr_synth=100 #nomber of molecule to generate
n_frag=10 #number of sequence fragments in each molecule
i_error=0.00 #insertion error rate
d_error=0.00 #deletion error rate
s_error=0.00 #substitution error rate

#----- parameters for sequencing -----#
gpu=false #use the gpu version of guppy basecaller, the cpu version will be used if false
nbr_read=400 #number of read
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
tag_size : $tag_size
spacer_path : $spacer_path
primers_path : $primers_path
nbr_synth : $nbr_synth
n_frag : $n_frag
i_error : $i_error
d_error : $d_error
s_error : $s_error
gpu : $gpu
nbr_read : $nbr_read
perfect_sequencing : $perfect

--results--
eof


#-----------------------------------------------------#
######### ===== Part 1: base sequences ====== #########
#-----------------------------------------------------#
start_time=$(date +"%s")

base_seq_path="$process_path/1_base_seq_file.fasta"

echo "Lancement du Workflow 4 : (taille $size_seq)"
echo "___Initialisation de la séquence de base___"

seq_gen_script="$project_dir/synthesis_simulation/sequence_generator/sequence_generator.py" #script for the sequence generation
if [ "$random_seq" = true ]
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
#-------------------------------------------------------------#
######### ===== Part 2: fragmentation with tag ====== #########
#-------------------------------------------------------------#
start_time=$(date +"%s")
echo "___Fragmentation des séquences___"

fragmentation_script="$project_dir/synthesis_simulation/fragmentation/fragmentation_with_tag_and_spacer.py" #script for the fragmentation
fragmentation_seq_path="$process_path/2_fragmented_seq_file.fasta"
python3 $fragmentation_script "$base_seq_path" "$fragmentation_seq_path" "$spacer_path" $frag_size $tag_size
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(date +"%s")
echo "fragmentation with tag : $(($end_time - $start_time)) ms" >> $time_file
#------------------------------------------------#
######### ===== Part 3: synthesis ====== #########
#------------------------------------------------#
start_time=$(date +"%s")
echo "___Synthèse des séquences___"

synthesis_script="$project_dir/synthesis_simulation/synthesis/molecule_synthesis.py" #script for the synthesis
synthesis_path="$process_path/3_synthesis_file.fasta"
python3 $synthesis_script -i "$fragmentation_seq_path" -o "$synthesis_path" -p $primers_path -n $nbr_synth --n_frag $n_frag --i_error $i_error --d_error $d_error --s_error $s_error
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

basecalling_path="$process_path/5_basecalling"
mkdir "$basecalling_path"

basecaller_path="$project_dir/sequencing_simulation/ont-guppy/bin/guppy_basecaller" #path of the basecaller to use

if [ "$gpu" = true ]
then
	echo "using GPU" 
	$basecaller_path -r --input_path "$sequencing_path" --save_path "$basecalling_path" -c dna_r9.4.1_450bps_hac.cfg --num_callers 1 --gpu_runners_per_device 1 --chunks_per_runner 256 --chunk_size 1024 -x "cuda:0"
else
	echo "using CPU"
	$basecaller_path -r --input_path "$sequencing_path" --save_path "$basecalling_path" -c dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 8 --num_callers 1
	
fi

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

reconstruction_path="$process_path/6_reconstruction"
mkdir "$reconstruction_path"

reconstruction_script="$project_dir/sequencing_simulation/spacer_sequencing/C++functions/src/reconstruction" #script for the reconstruction
fastq_file=(*.fastq)

$reconstruction_script "$basecalling_path/"$fastq_file "$reconstruction_path" "$spacer_path" $frag_size $tag_size

reconstruction_script_2="$project_dir/sequencing_simulation/spacer_sequencing/reconstruct_workflow_4_after_C.py" #script for the reconstruction
python3 $reconstruction_script_2 "$reconstruction_path" "$spacer_path" $frag_size $tag_size
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
	echo $line
    echo $line >> $summary
done

#-------------- Exit --------------#
echo "___Fin du processus \!___"

exit 0
