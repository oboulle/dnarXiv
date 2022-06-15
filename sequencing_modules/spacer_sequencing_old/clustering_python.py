import inspect
import os
import shutil
import sys
import networkx as nx
import consensus
import prog_dyn
import time
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0, parentdir+"/synthesis_simulation")

import utils.dna_file_reader as dfr
import utils.dna_numbering as dnbr


TRANSTAB = str.maketrans("ACGT", "TGCA")

# compute the reverse complement of a DNA sequence                                                                                                    
def reverseComplement(seq):
    global TRANSTAB
    return seq.translate(TRANSTAB)[::-1]

if len(sys.argv) != 6:
    print("usage : reconstruct.py input.fastq output_dir spacer.fasta fragment_size tag_size")
    sys.exit(1)

output_dir = sys.argv[2]
SPACER = dfr.read_fasta(sys.argv[3])["spacer"]
print ("Spacer:  ", SPACER)
FRAG_SIZE = int(sys.argv[4])
print ("Fragment size:", FRAG_SIZE)
TAG_SIZE = int(sys.argv[5])
print ("Tag size:  ", TAG_SIZE)


# read fastq sequences
SEQ = []
ff = open(sys.argv[1])
com = ff.readline()
while com != "":
    SEQ.append(ff.readline()[:-1])
    ff.readline()
    ff.readline()
    com = ff.readline()
ff.close()

print (len(SEQ),"sequences read")
start = time.time();

RC_SPACER = reverseComplement(SPACER)


WS = int(len(SPACER)/2) # size of the kmers

K_spacer = [] # kmers of the base spacer
K_rev_spacer = [] # kmers of reverse complement of the base spacer
for k in range(WS+1):
    kmer = SPACER[k:k+WS]
    K_spacer.append(kmer)
    kmer = RC_SPACER[k:k+WS]
    K_rev_spacer.append(kmer)
    
print ("Cut sequences into fragments")

FRAG = []
# cut sequences into fragments
for seq in SEQ:
    # locate spacer words
    L = []
    for j in range(len(seq)-WS+1):
        kmer = seq[j:j+WS]
        if kmer in K_spacer or kmer in K_rev_spacer:
            L.append([j,kmer])
    # group close spacer coordinates together
    COORD = []
    lastcoord = -100
    setcoord = []
    for couple in L:
        coord = couple[0]
        kmer = couple[1]
        if coord - lastcoord < int(len(SPACER)/2)+1:
            setcoord.append(couple)
        else:
            if len(setcoord)!= 0:
                COORD.append([])
                l = len(COORD)-1
                for c in setcoord:
                    COORD[l].append(c)
            setcoord = [couple]
        lastcoord = coord

    # determine spacer location
    for j in range (len(COORD)-1):
        end_spacer = COORD[j][-1] # spacer before the fragment
        start_next_spacer = COORD[j+1][0] # spacer after the fragment
        
        end_spacer_coord = end_spacer[0] + len(end_spacer[1]) # coord where the spacer before the fragment stop
        start_next_spacer_coord = start_next_spacer[0] # coord where the spacer after the fragment starts
        
        fragment = seq[end_spacer_coord:start_next_spacer_coord]
        
        if end_spacer[1] in K_rev_spacer and start_next_spacer[1] in K_rev_spacer: # the fragment has been reverted
            fragment = reverseComplement(fragment)
            
        FRAG.append(SPACER+fragment+SPACER) # re add the original spacer for the consensus
len_FRAG = len(FRAG);

print("cut sequences",time.time()-start, "s")

print (" -",len_FRAG,"fragments detected")

print ("Cluster fragments ")
# group similar fragments together
# --------------------------------
start = time.time();

# compute solid kmers
SKMER = {}
for fragment in FRAG:
    # count informative kmers (i.e. without spacers)
    for j in range(len(SPACER),len(fragment)-len(SPACER)-WS+1):
        kmer = fragment[j:j+WS]
        if kmer in SKMER:
            SKMER[kmer] += 1
        else:
            SKMER[kmer] = 1

# assign a set of kmers to each fragment
SOLID_THRESHOLD = 5
k=0
print("len dict",len(SKMER))
SET_FRAG = []
for fragment in FRAG:
    related_kmers = []
    for j in range(len(SPACER),len(fragment)-len(SPACER)-WS+1):
        kmer = fragment[j:j+WS]

        if SKMER[kmer] >= SOLID_THRESHOLD:
            k+=1
            related_kmers.append(kmer)
    SET_FRAG.append(related_kmers)
print("assign kmers",time.time()-start, "s")

print("kmer",k)
# construct fragment graph
G = nx.Graph()

#init a graph without edges
for i in range(len_FRAG):
    G.add_node(i)
k=0;
print(len(SET_FRAG))
start = time.time();

kmer_occurence = {} #each key is a kmer, the associated value is the list of index of fragments containing this kmer
for i in range(len_FRAG):
    for kmer in SET_FRAG[i]:
        if kmer in kmer_occurence:
            kmer_occurence[kmer].append(i)
        else:
            kmer_occurence[kmer] = [i]
            
for i in range(len_FRAG):
    print ("  ",int(i*100/len_FRAG),"%", end=chr(13))
    l = [0 for k in range(len_FRAG)] # l[k] = number of shared kmers with fragment k
    for kmer in SET_FRAG[i]:
        for num_frag in kmer_occurence[kmer]:
            l[num_frag]+=1

    for num_frag in range(i+1,len_FRAG):
        #if fragment i shares more than FRAG_SIZE/3 kmers with fragment k, an edge i-k is created in G
        #the connected components of G will be the clusters of fragments
        if l[num_frag] >= FRAG_SIZE/3:
            G.add_edge(i,num_frag)
            k+=1;

print("construct frag graph",time.time()-start, "s")
print("k : ",k)

#nx.write_gexf(G,"graph.gexf")  # debug
#G = nx.read_gexf("graph.gexf") # debug

k = 0
for cc in nx.connected_components(G):
    if len(cc) > 10:
        k += 1

print (" -",k,"clusters detected")
print ("Compute consensus fragments")

# compute consensus from connected components
# -------------------------------------------
CONSENSUS = []
k = 0
try:
    os.mkdir(output_dir+"/frag")  # create a fragment directory
except OSError:
    shutil.rmtree(output_dir+"/frag")
    os.mkdir(output_dir+"/frag")

for cc in nx.connected_components(G):
    if len(cc) > 10 :
        ff = open(output_dir+"/frag/f"+str(k)+".fasta","w")
        for n in cc:
            ff.write(">s"+str(k)+"_"+str(n)+"\n"+FRAG[int(n)]+"\n")
        ff.close()
        k += 1


