import sys
import random
import networkx as nx

def extractSolidKmer(PKMER,SOLID_THRESHOLD):
    SPKMER = {}
    SKMER  = {}
    maxpos = 0
    for pkmer in PKMER:
        w = int(pkmer.split(':')[2])
        if w >= SOLID_THRESHOLD:
            SPKMER[pkmer] = True
            kmer = pkmer.split(':')[0]
            pos = int(pkmer.split(':')[1])
            maxpos = max(maxpos,pos)
            if kmer not in SKMER:
                SKMER[kmer] = []
            SKMER[kmer].append(pkmer)
    return SKMER, SPKMER, maxpos

def searchPKmer(SKMER,kmer,ix):
    if kmer in SKMER:
        for pkmer in SKMER[kmer]:
            k = int(pkmer.split(':')[1])
            if abs(k-ix) < 10 :
                return pkmer
    return "none"

PATH  = {}
DONE  = {}
def searchPath(G,pkmer,start_pkmer):
    global SCORE
    global DONE
    global PATH
    if DONE[pkmer] == True :
        return
    DONE[pkmer] = True
    if pkmer == start_pkmer:
        l = len(start_pkmer.split(':')[0])
        PATH[start_pkmer][l] = [int(start_pkmer.split(':')[2]),start_pkmer.split(':')[0]] # score, sequence
        DONE[start_pkmer] = True
        return 
    pkmer_score = int(pkmer.split(':')[2])
    for pred_pkmer in G.predecessors(pkmer):
        edge_weight = G[pred_pkmer][pkmer]['weight']
        searchPath(G,pred_pkmer,start_pkmer)
        for l1 in PATH[pred_pkmer]:
            p = PATH[pred_pkmer][l1]
            score = p[0] + edge_weight*pkmer_score
            extend = pkmer.split(':')[0][-edge_weight:]
            path = p[1]+extend
            l2 = len(path)
            if l2 not in PATH[pkmer]:
                PATH[pkmer][l2] = [score,path]
            else:
                if score > PATH[pkmer][l2][0]:
                    PATH[pkmer][l2] = [score,path]                
    return 
    
def fragConsensus(FRAG,LEN_FRAG,SPACER):
    KMER_SIZE = int(len(SPACER)/2)
    # get kmers and their positions
    KMER = {}
    for i in range(len(FRAG)):
        frag = FRAG[i]
        for j in range(len(frag)-KMER_SIZE+1):
            kmer = frag[j:j+KMER_SIZE]
            if kmer in KMER:
                KMER[kmer].append(j)
            else:
                KMER[kmer] = [j]
    # separate identical kmers that are located at differents positions
    for kmer in KMER:
        KMER[kmer].sort()
        l1 = []
        l2 = [KMER[kmer][0]]
        for i in range(1,len(KMER[kmer])):
            if KMER[kmer][i]-KMER[kmer][i-1] < 10:
                l2.append(KMER[kmer][i])
            else:
                if len(l2) != 0:
                    l1.append(l2)
                l2 = []
        if len(l2) != 0:
            l1.append(l2)
        KMER[kmer] = l1
    # construct PKMER dictionary
    # a pkmer is a concatenation of kmer & position & weight
    # example: ACCTGAGT:23:10
    PKMER = {}
    for kmer in KMER:
        for i in range(len(KMER[kmer])):
            l1 = KMER[kmer][i]
            s = 0
            for k in l1:
                s += k
            s = int(s/len(l1))
            name = kmer+":"+str(s)+":"+str(len(l1))
            PKMER[name] = True
    # construct "good" pkmer overlap graph (KOG)
    T = int(len(FRAG)/3)
    BEST_PATH = {}
    while T > 6:
        T -= 1
        SKMER, SPKMER, maxpos = extractSolidKmer(PKMER,T)
        G = nx.DiGraph()
        for i in range(len(FRAG)):
            frag = FRAG[i]
            for j in range(len(frag)-KMER_SIZE+1):
                kmer1 = frag[j:j+KMER_SIZE]
                pkmer1 = searchPKmer(SKMER,kmer1,j)
                if pkmer1 in SPKMER:
                    for k in range(1,int(KMER_SIZE/2)):
                        kmer2 = frag[j+k:j+k+KMER_SIZE]
                        pkmer2 = searchPKmer(SKMER,kmer2,j+k)
                        if pkmer2 in SPKMER:
                            p1 = int(pkmer1.split(':')[1])
                            p2 = int(pkmer2.split(':')[1])
                            if p1 < p2:
                                if pkmer1 not in G:
                                    G.add_node(pkmer1)
                                if pkmer2 not in G:
                                    G.add_node(pkmer2)
                                G.add_edge(pkmer1,pkmer2,weight=k)


        start_node = searchPKmer(SKMER,SPACER[-KMER_SIZE:],0)
        if start_node not in G:
            continue
        end_node   = searchPKmer(SKMER,SPACER[:KMER_SIZE],maxpos)
        if end_node not in G:
            continue
        if nx.has_path(G,start_node,end_node) == False:
            continue

        findCorrectLength = False
        # search path in KOG
        for pkmer in G:
            DONE[pkmer] = False
            PATH[pkmer] = {}
        searchPath(G,end_node,start_node)
        for l in PATH[end_node]:
            if l == LEN_FRAG+len(SPACER):
                findCorrectLength = True
            if l not in BEST_PATH:
                BEST_PATH[l] = PATH[end_node][l]
            else:
                if PATH[end_node][l][0] > BEST_PATH[l][0]:
                    BEST_PATH[l] = PATH[end_node][l]
        if findCorrectLength == True:
            break
    #nx.write_gexf(G,"graph_consensus.gexf")
    bp = [1000,""]
    for l in BEST_PATH:
        #print(l)
        d = abs(LEN_FRAG+len(SPACER) - l)
        if d < bp[0]:
            bp = [d,BEST_PATH[l][1]]
    return bp[1][KMER_SIZE:-KMER_SIZE]

def main():
    SPACER = "AAAAAAAAAACCCCCCCCCC"
    SPACER = "AAAAAAAACCCCCCCC"

    LEN_FRAG = 200
    # read fragments
    LIST_FRAG = []
    ff = open(sys.argv[1])
    com = ff.readline()
    while com != "":
        seq = ff.readline()[:-1]
        LIST_FRAG.append(seq)
        com = ff.readline()
    ff.close()

    print(sys.argv[1],len(LIST_FRAG))

    seq = fragConsensus(LIST_FRAG,LEN_FRAG,SPACER)

    print (seq,len(seq))

if __name__ == "__main__":
    main()
