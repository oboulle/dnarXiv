#echo "-------------------------------------------"
#echo "|  Channel_decoder: LDPC BP DECODER       |"
#echo "-------------------------------------------"


####
#->INPUTS:
#	-$1:consensus path
#	-$2:expected length of decoded sequences without redondancy
#	-$3:output_path 
#	-$4:validity_result_path
#->OUTPUTS:
#	-output_path : decoded sequences output
#	-validity_result_path : results of validity checking
####


decoder_dir="$(dirname $0)" #directory containing this script
consensus_path="$1"
kInf=$2
output_path="$3"
validity_result_path="$4"

tmp_output_path="${output_path/./_tmp.}" #file the save the sequences before removing channel coding redundancy

H_Matrix="$decoder_dir"/matrix/matrix_K${kInf}N200.txt

printf "channel decoding...\n"
#decode the sequences
octave --silent "$decoder_dir"/BPDecoder/DNA_BP_LDPC_decoder.m "$consensus_path" "$tmp_output_path" "$H_Matrix"
if [ ! $? = 0 ]
then
	echo "error channel decoding"
	exit 1
fi
printf "\tcompleted !\n"	
	
#test the validity of decoded sequences
julia "$decoder_dir"/NB_LDPC_decoder_DNA.jl "$tmp_output_path" "$H_Matrix" "$validity_result_path"
	
#remove the redundancy
if [ "$(uname)" == "Darwin" ] #use gsed in macOS
then
	gsed -E '0~2s/(.{'$kInf'}).*/\1/' "$tmp_output_path" > "$output_path"
else
	sed -E '0~2s/(.{'$kInf'}).*/\1/' "$tmp_output_path" > "$output_path"
fi

#delete temp sequences file
rm "$tmp_output_path"
