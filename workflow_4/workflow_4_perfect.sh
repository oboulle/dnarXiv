#!/bin/bash

#--------------------------------------------------------------#
######### ===== Set the parameters in this part ====== #########
#--------------------------------------------------------------#

process_name="workflow_4_perfect_test/workflow_4" #name of the directory to save the generated files

#----- parameters for sequence generation -----#
random_seq=true #generate random sequences (true) or use an existing fasta file (false)
seq_path="workflow_4_fail/1_base_seq_file.fasta" #if random_seq is false, the sequences from this path are used (one .fasta file)
#else the sequences are generated with the following parameters
nbr_seq=1 #number of sequences
size_seq=3060 #size of the sequences
h_max=3 #maximum size for the homopolymeres

#----- parameters for fragmentation -----#
frag_size=200 #size of the fragments
tag_size=13 #size of the tag
spacer_path="spacer.fasta" #path to the spacer to use (.fasta file)

#----- parameters for synthesis -----#
primers_path="primers.fasta" #path to the primers to use (.fasta file)
nbr_synth=250 #nomber of molecule to generate
n_frag=10 #number of sequence fragments in each molecule
i_error=0.00 #insertion error rate
d_error=0.00 #deletion error rate
s_error=0.00 #substitution error rate

#----- parameters for sequencing -----#
nbr_read=250 #number of read
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
tag_size : $tag_size
spacer_path : $spacer_path
primers_path : $primers_path
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

echo "Lancement du Workflow 4 : (taille $size_seq)"
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
#-----------------------------------------------------#
######### ===== Part 4: convert to fastq ====== #########
#-----------------------------------------------------#
fasta_to_fastq_script="$project_dir/synthesis_simulation/utils/dna_file_reader.py" #script to convert fasta to fastq
fastq_file="$process_path/4_sequencing.fastq"
python3 $fasta_to_fastq_script "$synthesis_path" "$fastq_file"
#-----------------------------------------------------#
######### ===== Part 5: reconstruction ====== #########
#-----------------------------------------------------#
start_time=$(date +"%s")
echo "___Reconstruction___"
reconstruction_script="$project_dir/sequencing_simulation/spacer_sequencing/reconstruct_workflow_4.py" #script for the reconstruction
reconstruction_path="$process_path/5_reconstruction"
mkdir "$reconstruction_path"
python3 $reconstruction_script "$fastq_file" "$reconstruction_path" "$spacer_path" $frag_size $tag_size
if [ ! $? = 0 ]
then
	exit 1
fi
end_time=$(date +"%s")
echo "reconstruction : $(($end_time - $start_time)) s" >> $time_file
#------------------------------------------------------#
######### ===== Part 6: Result Analysis ====== #########
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
