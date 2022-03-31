#!/bin/bash

#script for xml metadata management in the containers

init_metadata_file() {
	#create a new metadata file for a new container
	#$1 metadata_file path
	#$2 container name (string)
	#$3 simulation (boolean)
	#$4 fragment length (int)
	#$5 assembly type (ordered/unordered)
	#$6 indexing (boolean)
	#$7 spacer (dna string)
	#$8 sequencing_technology (string)
	cat > "$1" <<- eof
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
	#$1 metadata_file path
	#$2 param name
	param_value=$(xmlstarlet sel -t -v "/container_metadata/$2" "$1")
	echo $param_value
}


set_container_param() {
	#set a main parameter from the container metadata
	#$1 metadata_file path
	#$2 param name
	#$3 param_value
	xmlstarlet edit -L\
  		--update "/container_metadata/$2" \
  		--value "$3" $1
}


add_document() {
	#add a new document to the container
	#$1 metadata_file path
	#$2 document index
	xmlstarlet edit -L -s "/container_metadata/documents" --type elem --name document \
		--insert  "/container_metadata/documents/document[not(@id)]" --type attr --name id --value $2 \
		$1
}

	
add_doc_param() {
	#add a new parameter to the document of the container
	#$1 metadata_file path
	#$2 document index
	#$3 param name
	#$4 param_value
	xmlstarlet edit -L -s "/container_metadata/documents/document[@id=$2]" --type elem --name $3 --value $4 $1
}


get_doc_param() {
	#get a parameter from a document
	#$1 metadata_file path
	#$2 document index
	#$3 param name
	#param_value=$(xmllint --xpath "string(//container_metadata/documents/document[@id=$2]/$3)" "$1")
	param_value=$(xmlstarlet sel -t -v //container_metadata/documents/document[@id=$2]/$3 "$1")
	
	echo $param_value
}

