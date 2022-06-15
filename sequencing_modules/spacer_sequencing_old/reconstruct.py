import inspect
import os
import shutil
import sys
import networkx as nx
import consensus
import prog_dyn
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0, parentdir+"/synthesis_simulation")

import utils.dna_file_reader as dfr

TRANSTAB = str.maketrans("ACGT", "TGCA")

# compute the reverse complement of a DNA sequence                                                                                                    
def reverseComplement(seq):
    global TRANSTAB
    return seq.translate(TRANSTAB)[::-1]


if len(sys.argv) != 6:
    print("usage : reconstruct.py input.fastq output.fasta spacer.fasta fragment_size overlap_min")
    sys.exit(1)

output_file = sys.argv[2]
SPACER = dfr.read_fasta(sys.argv[3])["spacer"]
print ("Spacer:  ", SPACER)
FRAG_SIZE = int(sys.argv[4])
print ("Fragment size:", FRAG_SIZE)
OVLP_MIN = int(sys.argv[5])
print ("Overlap min:  ", OVLP_MIN)


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

# put all sequence in forward direction
KF = {}
KR = {}
WS = int(len(SPACER)/2)

RC_SPACER = reverseComplement(SPACER)
for k in range(len(SPACER)-WS+1):
    kmer = SPACER[k:k+WS]
    KF[kmer] = k
    kmer = RC_SPACER[k:k+WS]
    KR[kmer] = k

for i in range(len(SEQ)):
    seq = SEQ[i]
    nbkf = 0
    nbkr = 0
    for j in range(len(seq)-WS+1):
        kmer = seq[j:j+WS]
        if kmer in KF:
            nbkf += 1
        if kmer in KR:
            nbkr += 1
    if nbkr > nbkf:
        SEQ[i] = reverseComplement(SEQ[i])


print ("Cut sequences into fragments")

FRAG = []
# cut sequences into fragments
for i in range(len(SEQ)):
    seq = SEQ[i]
    # locate spacer words
    L = []
    for j in range(len(seq)-WS+1):
        kmer = seq[j:j+WS]
        if kmer in KF:
            nbkf += 1
            L.append([j,kmer])

    # group close spacer coordinates together
    COORD = []
    lastcoord = -100
    setcoord = []
    for j in range(len(L)):
        coord = L[j][0]
        kmer = L[j][1]
        if coord - lastcoord < WS+1:
            setcoord.append(L[j])
        else:
            if len(setcoord)!= 0:
                COORD.append([])
                l = len(COORD)-1
                for c in setcoord:
                    COORD[l].append(c)
            setcoord = [L[j]]
        lastcoord = coord

    # determine spacer location
    SC = []
    for j in range (len(COORD)):
        if len(COORD[j]) > 2 :
            coord = COORD[j][0][0]
            kmer  = COORD[j][0][1]
            coord -= KF[kmer]
            SC.append(coord)

    # cut into fragments
    for j in range(len(SC)-1):
        l = SC[j+1]-SC[j] - len(SPACER)
        if l > FRAG_SIZE  - FRAG_SIZE/10 and l < FRAG_SIZE + FRAG_SIZE/10 :
            FRAG.append(seq[SC[j]:SC[j+1]+len(SPACER)])

print (" -",len(FRAG),"fragments detected")

print ("Cluster fragments ")
# group similar fragments together
# --------------------------------

# compute solid kmers
SKMER = {}
SOLID_THRESHOLD = 5
for i in range(len(FRAG)):
    frag = FRAG[i]
    # count informative kmers (i.e. without spacers)
    for j in range(len(SPACER),len(frag)-len(SPACER)-WS+1):
        kmer = frag[j:j+WS]
        if kmer in SKMER:
            SKMER[kmer] += 1
        else:
            SKMER[kmer] = 1
list_non_solid_kmer = []
for kmer in SKMER:
    if SKMER[kmer] < SOLID_THRESHOLD:
        list_non_solid_kmer.append(kmer)
for kmer in list_non_solid_kmer:
    del (SKMER[kmer])

# assign a set of kmers to each fragment
SET_FRAG = []
for i in range(len(FRAG)):
    x = set()
    frag = FRAG[i]
    for j in range(len(SPACER),len(frag)-len(SPACER)-WS+1):
        kmer = frag[j:j+WS]
        if kmer in SKMER:
            x.add(kmer)
    SET_FRAG.append(x)

# construct fragment graph
G = nx.Graph()
for i in range(len(FRAG)):
    G.add_node(i)
for i in range(len(FRAG)):
    print ("  ",int(i*100/len(FRAG)),"%",end=chr(13))
    for j in range(len(FRAG)):
        if i != j :
            d = len(SET_FRAG[i] & SET_FRAG[j])
            if d >= FRAG_SIZE/3:
                G.add_edge(i,j)

#nx.write_gexf(G,"graph.gexf")  # debug
#G = nx.read_gexf("graph.gexf") # debug

k = 0
for cc in nx.connected_components(G):
    if len(cc) > 10 :
        k += 1

print (" -",k,"clusters detected")
print ("Compute consensus fragments")

# compute consensus from connected components
# -------------------------------------------
CONSENSUS = []
k = 0
try:
    os.mkdir("frag")  # create a fragment directory
except OSError:
    shutil.rmtree("frag")
    os.mkdir("frag")

for cc in nx.connected_components(G):
    if len(cc) > 10 :
        ff = open("frag/f"+str(k)+".fasta","w")
        F = []
        for n in cc:
            F.append(FRAG[int(n)])
            ff.write(">s"+str(k)+"\n"+FRAG[int(n)]+"\n")
        seqc = consensus.fragConsensus(F,FRAG_SIZE,SPACER)
        CONSENSUS.append(seqc)
        print ("  ",k,len(seqc))
        ff.close()
        k += 1

# assemble fragments into one sequence
# ------------------------------------
print("Assemble fragments")
# overlap graph
OVLP = nx.DiGraph()
for i in range(len(CONSENSUS)):
    OVLP.add_node(i)
for i in range(len(CONSENSUS)):
    fi = CONSENSUS[i]
    for j in range(len(CONSENSUS)):
        fj = CONSENSUS[j]
        w = prog_dyn.compareFrag(fi[-OVLP_MIN:],fj[:OVLP_MIN],10)
        if w < 5 :
            OVLP.add_edge(i,j)
            break

#nx.write_gexf(OVLP,"ovlp_graph.gexf") # debug

# search first node, i.e.node without predecessors
start_node = -1
for n in OVLP:
    p  = [pred for pred in OVLP.predecessors(n)]
    if len(p) == 0:
        start_node = n
        break

# search last node, i.e.node without successors
end_node = -1
for n in OVLP:
    p  = [succ for succ in OVLP.successors(n)]
    if len(p) == 0:
        end_node = n
        break

for path in nx.all_simple_paths(OVLP, source=start_node, target=end_node):
    if len(path) == len(CONSENSUS):
        seq = CONSENSUS[path[0]]
        for i in range(1,len(path)):
            seq += CONSENSUS[path[i]][OVLP_MIN:]

ff = open(output_file,"w")
ff.write (">sequence built from "+sys.argv[1]+"\n"+seq+"\n")
ff.close()
print (" - sequence size =",len(seq))

