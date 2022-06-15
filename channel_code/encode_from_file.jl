#!/usr/local/bin/julia

#Apply LDPC channel encoding on a fasta file containing multiple sequences
#save the results in an output fasta file

################
#    Init()	   #
################

include("GF4Arithmetic.jl")

using MAT
using DelimitedFiles


###################
#    Functions    #
###################

function Dna2DecList(sequence)
	#convert a dna_sequence into a list of int

	sequence = replace(sequence, "A" => "0")
	sequence = replace(sequence, "C" => "1")
	sequence = replace(sequence, "G" => "2")
	sequence = replace(sequence, "T" => "3")
	dec_list = [parse(Int, char) for char in sequence]
	return dec_list
end


function DecList2Dna(dec_list)
	#convert a list of int into a dna sequence

	sequence = string(join(dec_list))
	sequence = replace(sequence, "0" => "A")
	sequence = replace(sequence, "1" => "C")
	sequence = replace(sequence, "2" => "G")
	sequence = replace(sequence, "3" => "T")
	return sequence
end


function apply_encoding(sequence)
	#apply the ldpc encoder on a sequence

	kInf = length(sequence)

	#get the abolute path of the matrix, independently from where the script has been called
	G_Matrix_path = string(@__DIR__) * "/matrix/G_matrix_K" * string(kInf) * "N200.txt"
	H_Matrix_path = string(@__DIR__) * "/matrix/matrix_K" * string(kInf) * "N200.txt"

	G = readdlm(G_Matrix_path)
	G=Int.(G)

	H = readdlm(H_Matrix_path)
	H = Int.(H) #use an integer matrix instead of float

	int_list=Dna2DecList(sequence)

	c=matrixGFMult(transpose(int_list),G)

	encoded_sequence=DecList2Dna(c)

	return encoded_sequence
end


################
#    Main()	   #
################

if(length(ARGS) != 2 )
	println("Usage: julia encode_from_file.jl source.fasta output.fasta")
	exit(1)
end

println("channel encoding...")

source_path=ARGS[1]
output_path=ARGS[2]

output = open(output_path, "w")

open(source_path) do source
	num_sequence = 0
	for i in enumerate(eachline(source))
		line = i[2]
		# skip the name of the sequences (starts with '>')
		if line[1] != '>'	
			encoded_line = apply_encoding(line)
			println(output, ">" * string(num_sequence))
			println(output, encoded_line)
			num_sequence += 1
		end
	end
end

close(output)

println("\tcompleted !")

