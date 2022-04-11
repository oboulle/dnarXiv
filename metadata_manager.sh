#!/bin/bash
set -u #exit and display error message if a variable is empty 

#script for xml metadata management in the containers
source "$(dirname ${BASH_SOURCE})"/log_manager.sh #load the log manager script

meta_file_name="metadata.xml" #name of the metadata file in containers

init_metadata_file() {
	#create a new metadata file for a new container
	#$1 container path
	#$2 container name (string)
	#$3 simulation (boolean)
	#$4 fragment length (int)
	#$5 assembly type (ordered/unordered)
	#$6 indexing (boolean)
	#$7 spacer (dna string)
	#$8 sequencing_technology (string)
	
	#save operation in logs
	log_metadata_call "$1" init_metadata_file "$2" "$3" "$4" "$5" "$6" "$7" "$8"
	
	cat > "$1"/$meta_file_name <<- eof
	<container_metadata>
	<!--Defined at the container creation-->
	  <name>"$2"</name>
	  <creation_date>$(date)</creation_date>
	  <simulation>$3</simulation>
	  <frag_length>$4</frag_length>
	  <assembly_type>$5</assembly_type>
	  <fragment_indexing>$6</fragment_indexing>
	  <spacer>$7</spacer>
	  <sequencing_technology>$8</sequencing_technology>
	<!--Can be modified when adding, deleting and storing documents-->
	  <number_of_documents>0</number_of_documents>
	  <editable>true</editable>
	  <documents><!--List of documents in the container-->
	  </documents>
	</container_metadata>
	eof
}


get_container_param() {
	#get a main parameter from the container me	tadata
	#$1 container path
	#$2 param name
	
	#save operation in logs
	log_metadata_call "$1" get_container_param "$2"
	
	xmlstarlet sel -t -v "/container_metadata/$2" "$1"/$meta_file_name
}


set_container_param() {
	#set a main parameter from the container metadata
	#$1 container path
	#$2 param name
	#$3 param_value
	
	#save operation in logs
	log_metadata_call "$1" set_container_param "$2" "$3"
	
	xmlstarlet edit -L\
  		--update "/container_metadata/$2" \
  		--value "$3" "$1"/$meta_file_name
}


add_document() {
	#add a new document to the container
	#$1 container path
	#$2 document index
	
	#save operation in logs
	log_metadata_call "$1" add_document "$2"
	
	#add a document then add the attribute 'index' for all documents wo attribute (should only be the added one)
	xmlstarlet edit -L -s "/container_metadata/documents" --type elem --name document \
		--insert  "/container_metadata/documents/document[not(@id)]" --type attr --name id --value "$2" \
		"$1"/$meta_file_name
	#add the creation date
	add_doc_param $1 $2 "creation_date" "$(date)"
	
	#----Update number of documents----#
	number_doc=$(get_container_param "$1" "number_of_documents")
	set_container_param "$1" "number_of_documents" $((number_doc+1))
}


del_document() {
	#delete a document from the container
	#$1 container path
	#$2 document index
	
	#save operation in logs
	log_metadata_call "$1" del_document "$2"
	
	xmlstarlet edit -L --delete "/container_metadata/documents/document[@id=$2]" $1/$meta_file_name
	
	#----Update number of documents----#
	number_doc=$(get_container_param "$1" "number_of_documents")
	set_container_param "$1" "number_of_documents" $((number_doc-1))
}

	
add_doc_param() {
	#add a new parameter to the document of the container
	#$1 container path
	#$2 document index
	#$3 param name
	#$4 param_value
	
	#save operation in logs
	log_metadata_call "$1" add_doc_param "$2" "$3" "$4"
	
	xmlstarlet edit -L -s "/container_metadata/documents/document[@id=$2]" --type elem --name "$3" --value "$4" "$1"/$meta_file_name
}


get_doc_param() {
	#get a parameter from a document
	#$1 container path
	#$2 document index
	#$3 param name
	
	#save operation in logs
	log_metadata_call "$1" get_doc_param "$2" "$3"
	
	param_value=$(xmlstarlet sel -t -v "/container_metadata/documents/document[@id=$2]/$3" "$1"/$meta_file_name)
	if [ "$param_value" != "" ]; then 
		echo "$param_value" 
	else
		echo "$3"
	fi
}


get_all_doc_param() {
	#get all parameters from a document
	#$1 container path
	#$2 document index
	
	#save operation in logs
	log_metadata_call "$1" get_all_doc_param "$2"
	
	xmlstarlet sel -T -t -m "/container_metadata/documents/document[@id=$2]/*" -v "concat(name(),'=',.,';')" "$1"/$meta_file_name
}




