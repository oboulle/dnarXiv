#!/usr/local/bin/julia

#---------------------------------------------------------------------------#
#	@Authors: Belaid Hamoum, Elsa Dupraz				    #
#	Project: DnarXiv, funded by CominLabs				    #
#---------------------------------------------------------------------------#
# 			DNA Data Storage SIMULATOR			    #
#      -Synthesis, Storage, preparation, nanopore sequencing, basecalling - #
#									    #
#---------------------------------------------------------------------------#

Bases=['A','C','G','T']

function insertB(seqTmp,rangeLenIns,lenIns)
	r3=rand()
	len=1;

	#length of the insertion
	for j=1:lenIns
		if (r3<=rangeLenIns[j])
			len=j;
			break;
		end
	end

	for j=1:len #insert after current pos
		base=rand(Bases) #pickup randomly a base to insert
		push!(seqTmp,base)
	end
end


function deleteB()
	len=1; #length of del is 1 (bursts are taken into account in another way)
	#Do nothing => deletion
end


function substB(seqTmp,rangeTransProb,currBase)
	r=rand()
	base=""
	function getSubstBase(indexB)
		for j=1:length(rangeTransProb[indexB])
			if (r<=rangeTransProb[indexB][j])
				base=split(subList[indexB,j],"2")[2] 
				break;
			end
		end
		return base
	end
	if(currBase=='A')
		base=getSubstBase(1)
	elseif(currBase=='C')
		base=getSubstBase(2)
	elseif(currBase=='G')
		base=getSubstBase(3)
	elseif(currBase=='T')
		base=getSubstBase(4)
	end
	push!(seqTmp,base)
end


function matchB(seqTmp,currBase)
	push!(seqTmp,currBase)
end


#####################
#	   Main()		#
#####################

using DelimitedFiles


if(length(ARGS) < 2 || length(ARGS) > 3)
	println("usage: script.jl input_path output_path k?")
	exit(1)
end

println("sequencing and basecalling simulation...")

input_path=ARGS[1]
output_path=ARGS[2]

k=6
if(length(ARGS) == 3)
	k=parse(Int64,ARGS[3])
end

simu_dir=string(@__DIR__) * "/probEdit/"

tmp=readdlm(simu_dir * "BegErrByPosAvg.txt");
BegProbEdit=Float64.(tmp[1]) 
tmp=readdlm(simu_dir * "BegInsByPosAvg.txt")
BegProbIns=Float64.(tmp[1])
tmp=readdlm(simu_dir * "BegDelByPosAvg.txt")
BegProbDel=Float64.(tmp[1])
tmp=readdlm(simu_dir * "BegMisAvg.txt")
BegProbSubst=Float64.(tmp[1])


MidProbEdit = Dict()
MidProbIns = Dict()
MidProbDel = Dict()
MidProbSubst = Dict()

for base in Bases
	tmp_prob=readdlm(simu_dir * base * "ErrAvg.txt");
	MidProbEdit[base] = Float64.(tmp_prob[1])
	
	tmp_prob=readdlm(simu_dir * base * "InsAvg.txt")
	MidProbIns[base] = Float64.(tmp_prob[1])
	
	tmp_prob=readdlm(simu_dir * base * "DelAvg.txt")
	MidProbDel[base] = Float64.(tmp_prob[1])
	
	tmp_prob=readdlm(simu_dir * base * "MisAvg.txt")
	MidProbSubst[base] = Float64.(tmp_prob[1])
end

subList=["A2C" "A2G" "A2T";
		 "C2A" "C2G" "C2T";
		 "G2A" "G2C" "G2T";
		 "T2A" "T2C" "T2G"]

Yi=["M" "D" "I" "S"]

transProb=[]
for i=1:length(subList[:,1])
	tmpProb=zeros(Float64,0)
	for j=1:length(subList[1,:])
		avgTmp=Float64.(readdlm(simu_dir * "$(subList[i,j])_Avg.txt")[1]);
		push!(tmpProb,avgTmp);
	end
	push!(transProb,tmpProb)
end

for i=1:length(subList[:,1])
	transProb[i]=transProb[i]/sum(transProb[i])
end

rangeTransProb=[] #2D table to allow easy mapping to subList

for i=1:length(subList[:,1])
	tmpRangeLine=[]
	for j=1:length(subList[1,:])
		if(j>1)
			push!(tmpRangeLine,tmpRangeLine[j-1]+transProb[i][j]);
		else
			push!(tmpRangeLine,transProb[i][j]);
		end
	end
	push!(rangeTransProb,tmpRangeLine)
end

#Probability to observe Yi on kmer i, knowing that Yi-1 occured
#index: 1-> prev Match, 2-> prev Del, 3-> prev Ins, 4-> prev Subst,
EditYiKmerprevYi=[]
mapKmerPrevYi=[]

#Yi-1= [Match[Yi=I,D,S,E,M],Del[Yi=I,D,S,E,M],Ins[Yi=I,D,S,E,M],Sub[Yi=I,D,S,E,M]]
for e1=1:4
	if(isfile(simu_dir * "KmerYi_prevYi$(Yi[e1])_RatesAvg.txt") )
		#kmer Ins Del Subst Err Match
		push!(EditYiKmerprevYi,readdlm(simu_dir * "KmerYi_prevYi$(Yi[e1])_RatesAvg.txt"));
		nKmer=length(EditYiKmerprevYi[e1][:,1]) #nbr of entries

		#compute avg probabilities (useful to use them when an unknown kmer is read)
		#add an extra index containing average probabilities that we can use when k-mer unknown
		avgIns=sum(EditYiKmerprevYi[e1][:,2])/nKmer
		avgDel=sum(EditYiKmerprevYi[e1][:,3])/nKmer
		avgSub=sum(EditYiKmerprevYi[e1][:,4])/nKmer
		avgErr=sum(EditYiKmerprevYi[e1][:,5])/nKmer
		avgMatch=sum(EditYiKmerprevYi[e1][:,6])/nKmer
		EditYiKmerprevYi[e1]=[EditYiKmerprevYi[e1] ; "" avgIns avgDel avgSub avgErr avgMatch];

		#create a dictionar (kmer,index) to reach the index (in KmerProb)containing probabilities (edit Rates table)
		push!(mapKmerPrevYi,Dict())
		nKmer=nKmer+1 #to take in account last added row (avgvalues)
		for i=1:nKmer
			mapKmerPrevYi[e1][EditYiKmerprevYi[e1][i,1]]=i
		end
	else
		push!(mapKmerPrevYi,Dict())
		push!(EditYiKmerprevYi,[]);
	end
end

#index: 1-> prev Match, 2-> prev Del, 3-> prev Ins, 4-> prev Subst,

kmerInsLenProbPrevYi=[]
mapKmerInsLenPrevYi=[]
#Yi-1= [Match[Yi=I,D,S,E,M],Del[Yi=I,D,S,E,M],Ins[Yi=I,D,S,E,M],Sub[Yi=I,D,S,E,M]]
for e1=1:4
	if(isfile(simu_dir * "KmerInsLen_prevYi$(Yi[e1])_RatesAvg2.txt") )
		kmerInsLenProbTmp=readdlm(simu_dir * "KmerInsLen_prevYi$(Yi[e1])_RatesAvg2.txt");

		#create a dictionar (kmer,index) to reach the index (in KmerProb)containing probabilities (edit Rates table)

		push!(mapKmerInsLenPrevYi,Dict())

		push!(kmerInsLenProbPrevYi,[]) #will contains prob of ins Length by column (each column => length) instead of string with ',' char
		nKmer=length(kmerInsLenProbTmp[:,1]) #nbr of entries
		for i=1:nKmer
			mapKmerInsLenPrevYi[e1][kmerInsLenProbTmp[i,1]]=i
			tmpLine=parse.(Float64,split(join(kmerInsLenProbTmp[i,2:end]),","))
			#do accumulation to divide ranges (1ins (0->0.20 ) , 2ins (0.20->0.30),...)
			for j=1:length(tmpLine)
				if(j>1)
						tmpLine[j]+=tmpLine[j-1]
				end
			end
			push!(kmerInsLenProbPrevYi[e1],tmpLine)
		end
		kmerInsLenProbTmp=Nothing; # remove tmp

		lenMax=length(kmerInsLenProbPrevYi[e1][1]); #number of column
		avgInsLen=[];
		#Avg value to get an insertion of length c
		for c=1:lenMax
			avg=0
			for k=1:nKmer
				avg+=kmerInsLenProbPrevYi[e1][k][c]
			end
			avg=avg/nKmer
			push!(avgInsLen,avg);
		end

		push!(kmerInsLenProbPrevYi[e1],Float64.(avgInsLen))
		#kmerInsLenProbPrevYi[e1]=[kmerInsLenProbPrevYi[e1];];
		#add an extra index to reach average probabilities that we can use when k-mer unknown
		mapKmerInsLenPrevYi[e1][""]=nKmer+1
	else
		push!(mapKmerInsLenPrevYi,Dict())
		push!(kmerInsLenProbPrevYi,[])
	end
end

kmerInsLenProbTmp=readdlm(simu_dir * "KmerInsLenRatesAvg2.txt");
#kmer probInsLen=1 probInsLen=2 probInsLen=3
#create a dictionar (kmer,index) to reach the index (in kmerInsLenProb)containing probabilities of insertion length
mapInsLenKmer=Dict()
kmerInsLenProb=[] #will contains prob of ins Length by column (each column => length) instead of string with ',' char
nKmer=length(kmerInsLenProbTmp[:,1]) #nbr of entries
for i=1:nKmer
	mapInsLenKmer[kmerInsLenProbTmp[i,1]]=i
	tmpLine=parse.(Float64,split(join(kmerInsLenProbTmp[i,2:end]),","))
	#do accumulation to divide ranges (1ins (0->0.20 ) , 2ins (0.20->0.30),...)
	for j=1:length(tmpLine)
		if(j>1)
			tmpLine[j]+=tmpLine[j-1]
		end
	end
	push!(kmerInsLenProb,tmpLine)
end

kmerInsLenProbTmp=Nothing # remove tmp

#probDel of current k-mer knowing that a delete occurs previously (Yi(kmeri) depend on Y-1)
DelProb_YiprevDel=readdlm(simu_dir * "KmerDel_delPrevRates.txt");
#kmer Occ Ins Del Subst Err Match
#create a dictionar (kmer,index) to reach the index (in KmerProb)containing probabilities (edit Rates table)
mapDelProb_YiprevDel=Dict()
nKmer=length(DelProb_YiprevDel[:,1]) #nbr of entries
for i=1:nKmer
	mapDelProb_YiprevDel[DelProb_YiprevDel[i,1]]=DelProb_YiprevDel[i,2]
end

tmp=readdlm(simu_dir * "EndErrByPosAvg.txt");
EndProbEdit=Float64.(tmp[1]) #majority vote Tab)
tmp=readdlm(simu_dir * "EndInsByPosAvg.txt")
EndProbIns=Float64.(tmp[1])
tmp=readdlm(simu_dir * "EndDelByPosAvg.txt")
EndProbDel=Float64.(tmp[1])
tmp=readdlm(simu_dir * "EndMisAvg.txt")
EndProbSubst=Float64.(tmp[1])

tmp=readdlm(simu_dir * "insLenBegRates.txt")[1,:]; #as a vector
pop!(tmp); #remove last blank
probInsLenBeg=Float64.(tmp) #majority vote Tab)

probInsLenMid = Dict()
for base in Bases
	tmp_prob = readdlm(simu_dir * base * "insLenMidRates.txt")[1,:];
	pop!(tmp_prob); #remove last blank
	probInsLenMid[base] = Float64.(tmp_prob) #majority vote Tab)
end

tmp=readdlm(simu_dir * "insLenEndRates.txt")[1,:];
pop!(tmp); #remove last blank
probInsLenEnd=Float64.(tmp) #majority vote Tab)

#divide space (0->1, depending on their rates) to pickup the type of edition (3 types)
#beg
sumBeg=BegProbIns+BegProbDel+BegProbSubst
rangeInsMaxBeg=BegProbIns/sumBeg;
rangeDelMaxBeg=rangeInsMaxBeg+(BegProbDel/sumBeg);
rangeSubstMaxBeg=rangeDelMaxBeg+(BegProbSubst/sumBeg);

#mid
sumMid = Dict()
rangeDelMaxMid = Dict()
rangeSubstMaxMid = Dict()
for base in Bases
	sumMid[base] = MidProbDel[base] + MidProbSubst[base]
	rangeDelMaxMid[base] = (MidProbDel[base]/sumMid[base]);
	rangeSubstMaxMid[base] = rangeDelMaxMid[base] + (MidProbSubst[base]/sumMid[base]);
end

#endrangeLenInsBeg
sumEnd=EndProbIns+EndProbDel+EndProbSubst
rangeInsMaxEnd=EndProbIns/sumEnd;
rangeDelMaxEnd=rangeInsMaxEnd+(EndProbDel/sumEnd);
rangeSubstMaxEnd=rangeDelMaxEnd+(EndProbSubst/sumEnd);

#divide  space to pickup the length of the edition
#Beg
lenInsBeg=length(probInsLenBeg);
rangeLenInsBeg=[]
for j=1:lenInsBeg	 
	if(j>1)
		push!(rangeLenInsBeg,rangeLenInsBeg[j-1]+probInsLenBeg[j])
	else
		push!(rangeLenInsBeg,probInsLenBeg[1])
	end
end

#mid

lenInsMid = Dict()
rangeLenInsMid = Dict()

for base in Bases
	lenInsMid[base] = length(probInsLenMid[base]);
	rangeLenInsMid[base]=[]
	for j=1:lenInsMid[base] 
		if(j>1)
			push!(rangeLenInsMid[base],rangeLenInsMid[base][j-1]+probInsLenMid[base][j])
		else
			push!(rangeLenInsMid[base],probInsLenMid[base][1])
		end
	end
end

lenInsEnd=length(probInsLenEnd);
rangeLenInsEnd=[]
for j=1:lenInsEnd	
	if(j>1)
		push!(rangeLenInsEnd,rangeLenInsEnd[j-1]+probInsLenEnd[j])
	else
		push!(rangeLenInsEnd,probInsLenEnd[1])
	end
end


function simulation(sequence)
	
	if !(last(sequence) in Bases)
		pop!(sequence) #remove last element (\n or '' ,...)
	end
	
	simSeq=[]

	index_seq=1	#used to go through the reference sequence
	kmer=""
	prevKmer=""   #used to know which kmer was before (to process ) NOT USED YET
	prevEdit=""   #used to know which edition was previously observed [1:M,2:D,3:I,4:S]

	if(k>1) #first base is affected by beg probabilities
		startKmer=1 #first Kmer pos
	else
		startKmer=2 #first Kmer pos
		#prevEdit="M" #As first Kmer will start at pos k first prevEdit will be at k pos
	end

	while index_seq <= length(sequence)

		currBase=sequence[index_seq]
		if(index_seq>1) #when 1<index_seq<k => use posBypos (length(kmer)<k =>posBypos done automatically )

			if(index_seq<=length(sequence)-1) 

				kmer=join(sequence[startKmer:index_seq]) 

				if(index_seq>=k ) #we have at least k base on the nanopore
					indexP=0

					if(get(mapKmerPrevYi[prevEdit],kmer,0)>0) #kmer entry exist
						indexP=mapKmerPrevYi[prevEdit][kmer] #mapKmer=["ACGA"]
					else #unknown kmer -> use average values
						indexP=mapKmerPrevYi[prevEdit][""]

						# (BACKtrack unknown k-mers)
						if(index_seq>=k ) #save them only one time
							#write(logUnknownKmers,kmer,"\n")
						end
					end

					#kmer INS	DEL	SUBST	ERR	MATCH

					probInsK=EditYiKmerprevYi[prevEdit][indexP,2]
					probDelK=EditYiKmerprevYi[prevEdit][indexP,3]
					probSubstK=EditYiKmerprevYi[prevEdit][indexP,4]
					probMatchK=EditYiKmerprevYi[prevEdit][indexP,6]
					probEditK=probDelK+probSubstK #Insertions can only appear after susbt or Match (otherwise prob not considered) 

					if( rand() <= probEditK) #edit the kmer

						rangeDelMaxK=probDelK/probEditK
						rangeSubstMaxK=rangeDelMaxK+(probSubstK/probEditK)

						r2=rand()
				
						if(r2<=rangeDelMaxK) #Deletion
										
							deleteB()
							prevEdit=2 #DEL						
						elseif(r2<=rangeSubstMaxK) #Substitution

							substB(simSeq,rangeTransProb,currBase)
							prevEdit=4 #Subs
							r3=rand()
							if (r3<=probInsK) #INS :In This version insertion appears after Match or Subst (=>extra char between it and next kmer)
								indexPLen=0
								if(get(mapKmerInsLenPrevYi[prevEdit],kmer,0)>0) #entry exist
									indexPLen=mapKmerInsLenPrevYi[prevEdit][kmer] 
								else  #entry doesn't exist -> use avg values (of all known kmers)
									indexPLen=mapKmerInsLenPrevYi[prevEdit][""]
								end

								rangeLenInsK=kmerInsLenProbPrevYi[prevEdit][indexPLen]
								lenInsK=length(rangeLenInsK)
								insertB(simSeq,rangeLenInsK,lenInsK)
								prevEdit=3 #INS
							end
						end
					else
						matchB(simSeq,currBase)
						prevEdit=1 #Match

						r3=rand()
						if (r3<=probInsK) #INS:#In This version insertion appears after Match or Subst 
							indexPLen=0
							if(get(mapKmerInsLenPrevYi[prevEdit],kmer,0)>0) #entry exist
								indexPLen=mapKmerInsLenPrevYi[prevEdit][kmer] #mapKmer=["ACGA"]
							else  #entry doesn't exist -> use avg values (of all known kmers)
								indexPLen=mapKmerInsLenPrevYi[prevEdit][""]
							end

							rangeLenInsK=kmerInsLenProbPrevYi[prevEdit][indexPLen]
							lenInsK=length(rangeLenInsK)
							insertB(simSeq,rangeLenInsK,lenInsK)
							prevEdit=3 #INS
						end
					end
				else #number of bases in the nanopore are less than k -> P(YI|XI) 
					if( rand() <= MidProbEdit[currBase]) #edit the nt

						r2=rand()

						if (r2<=rangeDelMaxMid[currBase])
							deleteB()
							prevEdit=2 #DEL
						elseif (r2<=rangeSubstMaxMid[currBase]) #substitution
							substB(simSeq,rangeTransProb,currBase)
							prevEdit=4 #SUBS

							r3=rand()
							if (r3<=MidProbIns[currBase]) #insertion
								insertB(simSeq,rangeLenInsMid[currBase],lenInsMid[currBase])
								prevEdit=3 #INS
							end
						end
					else
						matchB(simSeq,currBase)
						prevEdit=1 #Match

						r3=rand()
						if (r3<=MidProbIns[currBase]) #insertion
							insertB(simSeq,rangeLenInsMid[currBase],lenInsMid[currBase])
							prevEdit=3 #INS
						end
					end
				end
				if(index_seq>=k)
					prevKmer=kmer
					startKmer+=1
				end
			else #i==n last position (end part)

				if( rand() <= EndProbEdit) #edit the nt

					r2=rand()

					if (r2<=rangeInsMaxEnd) #insertion
						insertB(simSeq,rangeLenInsEnd,lenInsEnd)
						prevEdit=3 #INS
					elseif (r2<=rangeDelMaxEnd) #deletion
						deleteB()
						prevEdit=2 #DEL
					elseif (r2<=rangeSubstMaxEnd) #substitution
						substB(simSeq,rangeTransProb,currBase)
						prevEdit=4 #SUBS
					end
				else
					matchB(simSeq,currBase)
					prevEdit=1 #Match
				end
			end
		else #i==1 || i<k first position (Beg part)

			if( rand() <= BegProbEdit) #edit the nt

				r2=rand()

				if (r2<=rangeInsMaxBeg) #insertion
					insertB(simSeq,rangeLenInsBeg,lenInsBeg)
					prevEdit=3 #INS
				elseif (r2<=rangeDelMaxBeg) #deletion
					deleteB()
					prevEdit=2 #DEL
				elseif (r2<=rangeSubstMaxBeg) #substitution
					substB(simSeq,rangeTransProb,currBase)
					prevEdit=4 #SUBS
				end
			else
				matchB(simSeq,currBase)
				prevEdit=1 #Match  ************
			end
		end
		index_seq+=1;
	end
	return join(simSeq)
end

println("\tloaded !")

output = open(output_path, "w")

open(input_path) do source
	num_sequence = 0
	for i in enumerate(eachline(source))
		line = i[2]
		# skip the name of the sequences (starts with '>')
		if line[1] != '>'
			simSeq = simulation(line)
			println(output, "@simu_mol_" * string(num_sequence))
			println(output, simSeq)
			println(output,"+")
			println(output,"---")
			num_sequence += 1
		end
	end
end

close(output)

println("\tcompleted !")
