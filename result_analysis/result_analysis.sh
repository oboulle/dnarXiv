

container_path="${1}"
document_index="${2}"

if [ "$#" != 2 ]
then
	echo "Usage : resulat_analysis container_path document_index"
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

result_analysis_script=$project_dir/fasta36/bin/ggsearch36
container_doc_path=$container_path/$document_index
$result_analysis_script -3 -O $container_doc_path/result.txt $container_doc_path/source.fasta $container_doc_path/reconstructed_source.fasta
