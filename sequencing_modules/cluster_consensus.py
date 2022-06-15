#!/usr/bin/python3

import os
import sys
import networkx as nx
import time

# used to get the consensus fragments from a directory of fragments clusters

class Pkmer:
    """
    class to save a kmer, its average position in a cluster and its weight
    """
    def __init__(self, kmer, position, weight):
        self.kmer = kmer
        self.position = position
        self.weight = weight


def extractSolidKmer(pkmer_list, SOLID_THRESHOLD):
    """
    pkmer = object : kmer,position,weight
    select the pkmers that have a weight > threshold    
    """
    solid_kmers  = {}
    solid_pkmers = {}
    maxpos = 0
    for pkmer in pkmer_list:
        if pkmer.weight >= SOLID_THRESHOLD:
            solid_pkmers[pkmer] = True
            
            if pkmer.kmer not in solid_kmers:
                solid_kmers[pkmer.kmer] = []
            solid_kmers[pkmer.kmer].append(pkmer)
            
            maxpos = max(maxpos,pkmer.position)
    
    return solid_kmers, solid_pkmers, maxpos


def searchPKmer(solid_kmers, kmer, frag_index):
    """
    selected_kmers : dict of (kmer:[pkmers])
    frag_index : position of the kmer in the fragment
    returns a pkmer of a kmer close to frag_index in the fragment
    """
    if kmer in solid_kmers:
        for pkmer in solid_kmers[kmer]:
            if abs(pkmer.position-frag_index) < 10 :
                return pkmer
    return "none"


PATH  = {} # dict of (pkmers, value = dict of (path_lengths, value=[score, sequence]) )
visited_pkmer  = {} # dict of pkmers, value = true if visited, false otherwise

def searchPath(G, pkmer, start_pkmer):
    global visited_pkmer
    global PATH
    
    if visited_pkmer[pkmer]:
        return
    visited_pkmer[pkmer] = True
    
    if pkmer == start_pkmer:
        len_kmer = len(start_pkmer.kmer)
        PATH[start_pkmer][len_kmer] = [start_pkmer.weight, start_pkmer.kmer] # score, sequence
        return
    
    for pred_pkmer in G.predecessors(pkmer):
        edge_weight = G[pred_pkmer][pkmer]['weight']
        searchPath(G, pred_pkmer, start_pkmer)
        
        for pred_path_length in PATH[pred_pkmer]:
            [pred_score, pred_path] = PATH[pred_pkmer][pred_path_length]
            
            score = pred_score + edge_weight*pkmer.weight
            extend = pkmer.kmer[-edge_weight:]
            path = pred_path+extend
            new_path_length = len(path)
            
            # if no paths of this length, or previous path of this length has lower score
            if new_path_length not in PATH[pkmer] or score > PATH[pkmer][new_path_length][0]:
                PATH[pkmer][new_path_length] = [score, path]      
    return 


def fragConsensus(frag_list, frag_length):
    """
    frag_length : original size of the fragment
    returns the consensus sequence from a list of fragments 
    """
    KMER_SIZE = 8 # arbitrary constant
    
    # small part to surround the fragments, help to identify the start and end, 
    # contains fictional base B that cannot exist in any fragments and should not be in the resulting consensus
    surrounding = KMER_SIZE*"B" 
    
    # get kmers and their positions
    kmers_positions = {} # dict of (kmer:[indexes list])
    for frag in frag_list:
        surr_frag = surrounding+frag+surrounding
        for i in range(len(surr_frag)-KMER_SIZE+1):
            kmer = surr_frag[i:i+KMER_SIZE]
            if kmer in kmers_positions:
                kmers_positions[kmer].append(i)
            else:
                kmers_positions[kmer] = [i]
    
    # separate identical kmers that are located at different positions
    for kmer in kmers_positions:
        kmers_positions[kmer].sort()
        clusters_list = [] #list of groups of close positions for the same kmer
        cluster = [kmers_positions[kmer][0]] # group of close positions for the same kmer
        for i in range(1,len(kmers_positions[kmer])):
            if kmers_positions[kmer][i]-kmers_positions[kmer][i-1] < 10:
                cluster.append(kmers_positions[kmer][i])
            else:
                if len(cluster) != 0:
                    clusters_list.append(cluster)
                cluster = []
        if len(cluster) != 0:
            clusters_list.append(cluster)
        kmers_positions[kmer] = clusters_list
    
    # construct a list of pkmers
    # a pkmer is an object containing a kmer & his average_position & weight
    # example: ACCTGAGT:23:10
    pkmer_list = []
    for kmer in kmers_positions:
        for cluster in kmers_positions[kmer]:
            s = sum(cluster)
            average_position = int(s/len(cluster))
            pkmer = Pkmer(kmer, average_position, len(cluster))
            pkmer_list.append(pkmer)
    
    # construct "good" pkmer overlap graph (KOG)
    THRESHOLD = int(len(frag_list)/3)
    BEST_PATHS = {} # dict of bests paths lengths, value=[score, sequence]
    while THRESHOLD > 6:
        THRESHOLD -= 1
        solid_kmers, solid_pkmers, maxpos = extractSolidKmer(pkmer_list, THRESHOLD)
        G = nx.DiGraph()
        for frag in frag_list:
            surr_frag = surrounding+frag+surrounding
            for i in range(len(surr_frag)-KMER_SIZE+1):
                kmer1 = surr_frag[i:i+KMER_SIZE]
                pkmer1 = searchPKmer(solid_kmers, kmer1, i)
                if pkmer1 in solid_pkmers:
                    for j in range(1,int(KMER_SIZE/2)):
                        kmer2 = surr_frag[i+j:i+j+KMER_SIZE]
                        pkmer2 = searchPKmer(solid_kmers, kmer2, i+j)
                        if pkmer2 in solid_pkmers:
                            if pkmer1.position < pkmer2.position:
                                if pkmer1 not in G:
                                    G.add_node(pkmer1)
                                if pkmer2 not in G:
                                    G.add_node(pkmer2)
                                G.add_edge(pkmer1, pkmer2, weight=j) # weight = number of shifts between the 2 kmers

        start_node = searchPKmer(solid_kmers, surrounding, 0)
        if start_node not in G:
            continue
        end_node = searchPKmer(solid_kmers, surrounding, maxpos)
        if end_node not in G:
            continue
        if not nx.has_path(G, start_node, end_node):
            continue

        # search path in kmer overlap graph
        for pkmer in G:
            visited_pkmer[pkmer] = False
            PATH[pkmer] = {}
        searchPath(G, end_node, start_node)
        
        findEspectedLength = False

        for path_length in PATH[end_node]:
            if path_length == frag_length+2*KMER_SIZE:
                findEspectedLength = True
            if path_length not in BEST_PATHS:
                BEST_PATHS[path_length] = PATH[end_node][path_length]
            else:
                if PATH[end_node][path_length][0] > BEST_PATHS[path_length][0]:
                    BEST_PATHS[path_length] = PATH[end_node][path_length]
                    
        if findEspectedLength == True:
            break
    
    #nx.write_gexf(G,"graph_consensus.gexf")
    best_path = [1000,""] #difference with expected length, consensus
    for path_length in BEST_PATHS:
        diff = abs(frag_length+2*KMER_SIZE - path_length)
        if diff < best_path[0]:
            best_path = [diff, BEST_PATHS[path_length][1]]
    return best_path[1][KMER_SIZE:-KMER_SIZE]


def consensusOnFragDir(frag_dir_path, output_path, frag_length):
    """
    apply the consensus algorithm on every cluster files from a directory
    """
    consensus_file = open(output_path,"w")
    computed_consensus = 1
    for cluster_filename in os.listdir(frag_dir_path):
        #display percentage progression
        i=int(20*computed_consensus/nbr_cluster)
        sys.stdout.write('\r')
        sys.stdout.write("[%-20s] %d%%" % ('='*i, 5*i))
        sys.stdout.flush()
        
        # read fragments
        frag_list = []
        cluster_file = open(os.path.join(frag_dir_path, cluster_filename))
        line = cluster_file.readline()
        while line != "":
            frag = cluster_file.readline()[:-1]
            frag_list.append(frag)
            line = cluster_file.readline()
        cluster_file.close()
        
        consensus_sequence = fragConsensus(frag_list, frag_length)
        if len(consensus_sequence) >= 2*frag_length/3:
            
            consensus_file.write(">consensus_of_"+cluster_filename+"\n")
            consensus_file.write(consensus_sequence+"\n")

        computed_consensus += 1
    consensus_file.close()


if __name__ == "__main__":
    
    if len(sys.argv) != 4:
        print("usage : consensus.py frag_dir output_path fragment_length")
        exit(1)
    
    print("consensus...")
    frag_dir_path = sys.argv[1]
    output_path = sys.argv[2]
    frag_length = int(sys.argv[3])

    nbr_cluster = len(os.listdir(frag_dir_path))
    if nbr_cluster == 0:
        print("error : no cluster file found")
        exit(1)
    #start = time.time()
    
    consensusOnFragDir(frag_dir_path, output_path, frag_length)
    
    #print("\n",time.time() - start,"seconds")
            
    print("\n\tcompleted !")

