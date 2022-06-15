//============================================================================
// Name        : Reconstruction.cpp
// Author      : Olivier
//============================================================================

//g++ -std=c++11 clustering.cpp -o clustering

#include <iostream>

#include <map>
#include <algorithm>
#include <vector>
#include <string>
#include "graph.h"
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <set>
#include <chrono>

using namespace std;

string read_spacer(string spacer_path);
vector<string> read_sequences(string sequences_path);
void clustering(vector<string> sequences, string spacer, int fragment_size, string output_dir);

int main(int argc,char* argv[]) {
	//args = input.fastq output_dir spacer.fasta fragment_size
	if(argc!=5){
		printf("usage : reconstruction.cpp input.fastq output_dir spacer.fasta fragment_size");
		return 1;
	}

	string input_path = argv[1]; //fastq file of the input sequences
	string output_dir = argv[2]; //directory to save the results
	string spacer_path = argv[3]; //path of the fasta spacer sequence

	string spacer = read_spacer(spacer_path);
	vector<string> sequences = read_sequences(input_path);
	int fragment_size = stoi(argv[4]);

	printf("Spacer:  %s\n", spacer.c_str());
	printf("Fragment size:  %d\n", fragment_size);

	vector<string> SEQ_1 = vector<string>(sequences.begin(), sequences.begin() + sequences.size()/2); // @suppress("Symbol is not resolved")
	vector<string> SEQ_2 = vector<string>(sequences.begin() + sequences.size()/2, sequences.end()); // @suppress("Symbol is not resolved")
	//apply a clustering to the first half of the sequences and another independent clustering to the second half
	clustering(SEQ_1, spacer, fragment_size, output_dir+"/frag1");
	clustering(SEQ_2, spacer, fragment_size, output_dir+"/frag2");
	printf("done !");

	return 0;
}

string read_spacer(string spacer_path) {
	//read the spacer sequence from a fasta file
	ifstream file(spacer_path);
	string line;
	getline(file, line);
	getline(file, line);
	file.close();
	return line;
}

vector<string> read_sequences(string sequences_path) {
	//read all the sequences from a fastq file
	vector<string> sequences;
	ifstream file(sequences_path);
	string line;
	while (getline(file, line)) {
		getline(file, line);
		sequences.push_back(line.c_str());
		getline(file, line);
		getline(file, line);
	}
	file.close();
	return sequences;
}

string reverse_complement(string sequence) {
	//return the reverse complement of a sequence
	reverse(sequence.begin(), sequence.end());
	for (size_t i = 0; i < sequence.length(); ++i){
		switch (sequence[i]){
	        case 'A':
	        	sequence[i] = 'T';
	            break;
	        case 'C':
	        	sequence[i] = 'G';
	            break;
	        case 'G':
	        	sequence[i] = 'C';
	            break;
	        case 'T':
	        	sequence[i] = 'A';
	            break;
	    }
	}
	return sequence;
}

void clustering(vector<string> sequences, string spacer, int fragment_size, string output_dir) {
	//indexation table clustering

	auto start = chrono::system_clock::now();

	int WS = int(spacer.size()/2); //size of the kmers

	set<string> K_spacer; //set of kmers from the spacer
	set<string> K_rev_spacer; //set of kmers from the reverse complement spacer
	string RC_spacer = reverse_complement(spacer);
	for(int k=0; k<WS+1; k++) {
		string kmer = spacer.substr(k, WS);
		K_spacer.insert(kmer);
		kmer = RC_spacer.substr(k, WS);
		K_rev_spacer.insert(kmer);
	}
	vector<string> fragments_vector;
	for(const string& sequence : sequences) {
		// locate spacer words
		vector<tuple<int, string>> spacer_kmers_coordinates; //contains pairs coord;kmer of spacers kmers in the sequence
		int end_range = sequence.size()-WS+1;
		for(int j=0; j<end_range; j++) {
			string kmer = sequence.substr(j, WS);
			if(K_spacer.find(kmer) != K_spacer.end()
					or K_rev_spacer.find(kmer) != K_rev_spacer.end()) {
				spacer_kmers_coordinates.push_back(make_tuple(j, kmer)); // @suppress("Invalid arguments")
			}
		}
		// group close spacer coordinates together
		vector<vector<tuple<int, string>>> COORD;
		int last_coord = -100;
		vector<tuple<int, string>> set_coord;
		for(const tuple<int, string>& couple : spacer_kmers_coordinates) {
			int coord = get<0>(couple);
			string kmer = get<1>(couple);
			if(coord - last_coord < int(spacer.size()/2)+1) {
				set_coord.push_back(couple);
			} else {
				if(set_coord.size() != 0) {
					COORD.push_back(vector<tuple<int, string>>());
					int l = COORD.size()-1;
					for(const tuple<int, string>& c : set_coord) {
						COORD[l].push_back(c);
					}
				}
				set_coord = vector<tuple<int, string>> {couple};
			}
			last_coord = coord;
		}
		//determine spacer location
		for(int j=0; j<COORD.size()-1; j++) {
			tuple<int, string> end_spacer = COORD[j].back(); //spacer before the fragment
			tuple<int, string> start_next_spacer = COORD[j+1][0]; //spacer after the fragment

			int end_spacer_coord = get<0>(end_spacer) + get<1>(end_spacer).size(); //coord where the spacer before the fragment stop
			int start_next_spacer_coord = get<0>(start_next_spacer); //coord where the spacer after the fragment starts

			string fragment = sequence.substr(end_spacer_coord,start_next_spacer_coord-end_spacer_coord);
			if(K_rev_spacer.find(get<1>(end_spacer)) != K_rev_spacer.end()
									and K_rev_spacer.find(get<1>(start_next_spacer)) != K_rev_spacer.end()) {
				//the fragment has been reverted
				fragment = reverse_complement(fragment);
			}
			fragments_vector.push_back(spacer+fragment+spacer); //re add the original spacer to the fragment for the consensus
		}
	}
	auto end = chrono::system_clock::now();

	chrono::duration<double> elapsed_seconds = end-start; // @suppress("Invalid arguments")
	time_t end_time = chrono::system_clock::to_time_t(end);

	cout << "cut sequences: " << elapsed_seconds.count() << "s\n";

	//ASSIGN SET KMER
	//assign a set of kmer to each fragment
	start = chrono::system_clock::now();
	cout << "assign_set_kmer " << fragments_vector.size() << endl;

	int len_spacer = spacer.size();
	map<string, int> list_kmers; //map of pairs kmers;nbr_occurrence in all the fragments

	for(const string& fragment : fragments_vector) {
		int end_range = fragment.size()-len_spacer-WS+1;
		string kmer;
		for(int j=len_spacer; j<end_range; j++) {
			kmer = fragment.substr(j, WS);
			list_kmers[kmer]++;
		}
	}

	int SOLID_THRESHOLD = 5; //minimum number of occurrence for a kmer in list_kmers

	vector<vector<string>> set_kmers_per_fragment; //list of kmers set, index = fragment index in the list of fragments
	for(const string& fragment : fragments_vector) {
		vector<string> related_kmers; //kmers of the fragment
		int end_range = fragment.size()-len_spacer-WS+1;
		string kmer;
		for(int j=len_spacer; j<end_range; j++) {
			kmer = fragment.substr(j, WS);
			if(list_kmers[kmer] >= SOLID_THRESHOLD){
				related_kmers.push_back(kmer);
			}
		}
		set_kmers_per_fragment.push_back(related_kmers);
	}

	end = chrono::system_clock::now();

	elapsed_seconds = end-start;
	end_time = chrono::system_clock::to_time_t(end);

	cout << "assign kmers: " << elapsed_seconds.count() << "s\n";


	//CONSTRUCT FRAGMENT GRAPH / clustering
	start = chrono::system_clock::now();

	int len_frag = fragments_vector.size();
	cout << "construct_fragment_graph " << set_kmers_per_fragment.size() << endl;
	Graph G(len_frag); //init a graph, node i = fragment of index i
	map<string, vector<int>> kmer_occurence; //map of pairs kmer;(list of index for fragments containing this kmer)

	for(int i=0; i<len_frag; i++) {
		for(const string& kmer : set_kmers_per_fragment[i]){
			kmer_occurence[kmer].push_back(i);
		}
	}

	for(int i=0; i<len_frag; i++) {
		vector<int> l(len_frag, 0); //vector where l[i] = number of shared kmers with fragment i
		for(const string& kmer : set_kmers_per_fragment[i]){
			for(const int& num_frag : kmer_occurence[kmer]){
				l[num_frag]++;
			}
		}
		for(int num_frag=i+1; num_frag<len_frag; num_frag++){
			if(l[num_frag]>=fragment_size/3){
				//an edge is created between 2 fragments if enough kmers are shared (fragment_size/3 kmers)
				G.add_edge(i, num_frag);
			}
		}
	}

	end = chrono::system_clock::now();

	elapsed_seconds = end-start;
	end_time = chrono::system_clock::to_time_t(end);

	cout << "construct frag graph: " << elapsed_seconds.count() << "s\n";

	//PRINT CONNECTED COMPONENTS
	//save the list of clusters fragments
	cout << "print_connected_components" << endl;
	int k = 0;
    cout << "mkdir "<< output_dir+" : " << mkdir((output_dir).c_str(), 0777) << endl;
    vector<vector<int>> connected_components = G.connected_components(G.v);
	cout << "connected components " << connected_components.size() << endl;
    for(const vector<int>& cluster : connected_components) {
    	//cout << cc.size() << endl;
    	if(cluster.size() > 10) {
    		ofstream outfile(output_dir+"/f"+to_string(k)+".fasta");
    		for(const int& node : cluster) {
    			//cout << node << ".";
    			outfile << ">s"+to_string(k)+"_"+to_string(node)+"\n"+fragments_vector[node] << endl;
    		}
    		//cout << endl;
    		outfile.close();
    		k++;
    	}
    }
	cout << k << " clusters" << endl;
}

int maint(int argc,char* argv[]) {
	/*add_kmer_to_dict("test");
	add_kmer_to_dict("test2");
	add_kmer_to_dict("test");
	add_kmer_to_dict("test3");*/

	//args = input.fastq output_dir spacer.fasta fragment_size tag_size
	if(argc!=6){
		printf("usage : reconstruction.cpp input.fastq output_dir spacer.fasta fragment_size tag_size");
		return 1;
	}

	string input_path = argv[1];
	string output_dir_path = argv[2];
	string spacer_path = argv[3];

	mkdir((output_dir_path+"/frag").c_str(), 0777);
	cout << "done" << endl;
	return 0;
}
