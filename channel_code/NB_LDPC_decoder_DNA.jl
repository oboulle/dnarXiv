#!/usr/local/bin/julia

#test the validity of a sequences file
#return a result file containing 0 for valid sequences and 1 for non valid sequences


################
#	Init()	   #
################

include("GF4Arithmetic.jl")

using MAT
using DelimitedFiles



###################
#    Functions    #
###################

#**connect v-nodes to c-nodes (each c-node will have indexes of all its v-nodes)**#
function VToCnodes(H,cnodes)

	n=length(H[:,1]) #number of line or equations or c-nodes
	
	for i=1:n
		push!(cnodes,[]) #new equation or connexion 
		m=length(H[i,:])
		for j=1:m
			if(H[i,j]!=0) #so the v-nodes is include on the parity check equation
				push!(cnodes[i],j) 
			end
		end
	end	
end


#check if code-word valid
function valid2(y,H,cnodes)
		
	n=length(cnodes)
	for i=1:n

		m=length(cnodes[i])
		tmp=0 
		for j=1:m 
			tmpM=multGF(y[cnodes[i][j]],H[i,cnodes[i][j]]) #
			tmp=addGF(tmp,tmpM); #	
		end
		if(tmp!=0)
			#code-word not valid
			return false
		end
	end
	return true
end

		
function Dna2DecList(sequence)
	#convert a dna_sequence into a list of int

	sequence = replace(sequence, "A" => "0")
	sequence = replace(sequence, "C" => "1")
	sequence = replace(sequence, "G" => "2")
	sequence = replace(sequence, "T" => "3")
	dec_list = [parse(Int, char) for char in sequence]
	return dec_list
end


################
#	Main()	   #
################

if(length(ARGS) != 3 )
	println("Usage: julia NB_LDPC_decoder_DNA.jl sequence_path H_matrix_path result_path")
	exit(1)
end

println("validity checking...")

sequence_path=ARGS[1]
matName=ARGS[2]
result_path=ARGS[3]

ext=split(matName,'.')[2]

H=[] #parity check matrix

if(ext=="mat")
	file = matopen(matName)
	H = Int.(read(file, "H")) 
else
	H = readdlm(matName)
	H = Int.(H) 
end

m=length(H[1,:]) #columns
n=length(H[:,1]) #rows


cnodes=[]
vnodes=[]

VToCnodes(H,cnodes) #init connexions v-nodes <-> cnodes

result = open(result_path, "w")

#read the sequences and check the validity
open(sequence_path) do source
	num_sequence = 0
	for i in enumerate(eachline(source))
		line = i[2]
		# skip the name of the sequences (starts with '>')
		if line[1] != '>'

			seq_len=length(line)
				
			if(seq_len!=m) 
				println(result, "1") #code-word not valid
			else
				
				vnodes=Dna2DecList(line)
				if(valid2(vnodes,H,cnodes))
					println(result, "0")	#code-word valid
				else
					println(result, "1")	#code-word not valid
				end
			end
			num_sequence += 1
		end
	end
end

println("\tcompleted !")
