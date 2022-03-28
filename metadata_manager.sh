#!/bin/bash



init_metadata_file() {
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
	  <simulation>$3</simulation>
	  <frag_length>$4</frag_length>
	  <assembly_type>$5</assembly_type>
	  <fragment_indexing>$6</fragment_indexing>
	  <spacer>$7</spacer>
	  <sequencing_technology>$8</sequencing_technology>
	
	<!--Can be modified when adding, deleting and storing documents-->
	  <number_of_documents>0</number_of_documents>
	  <editable>true</editable>
	  
	</container_metadata>
	eof
}


get_container_param() {
	#get a main parameter from the container me	tadata
	#$1 param name
	#$2 metadata_file path
	param_value=$(xmllint --xpath "string(//container_metadata/$1)" "$2")
	echo $param_value
}

set_container_param() {
	#set a main parameter from the container metadata
	#$1 param name
	#$2 param_value
	#$3 metadata_file path
	echo -e "cd //container_metadata/$1\nset $2\nsave"|xmllint --shell $3 >/dev/null
	
	#>/dev/null to hide the xmllint shell commands
}

get_doc_param() {
	#get a parameter from a document
	#$1 document index
	#$2 param name
	#$3 metadata_file path
	param_value=$(xmllint --xpath "string(//container_metadata/documents/document[@id=$1]/$2)" "$3")
	echo $param_value
}

init_metadata_file "meta_test.xml" "container_name" true 100 "ordered" false "AACGGTCGTC" "nanopore"

