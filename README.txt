Dnarxiv commands

dna_create
	dna_create [-sim] [-fl int] [-sp string] Cname
	create an empty container

dna_add
	dna_add Dname Cname
	add a document in the container

dna_del
	dna_del Cname DI
	delete a document from the container

dna_list
	dna_list Cname
	display the list of documents in the container

dna_store
	dna_store [-n_synth int] [-i_error int] [-d_error int] [-s_error int] Cname
	synthetise the documents of the container into dna molecules

dna_read
	dna_read Cname DI Dname
	sequence a document of the container

dna_workflow
	dna_workflow
	complete cycle test with a small document
