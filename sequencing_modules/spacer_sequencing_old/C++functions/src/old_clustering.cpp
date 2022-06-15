//============================================================================
// Name        : Reconstruction.cpp
// Author      : Olivier
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

//g++ -std=c++11 old_clustering.cpp -o old_clustering

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
void reconstruction(vector<string> SEQ, string SPACER, int FRAG_size, string output_dir);


map<string, string*> kmer_dict;

int main(int argc,char* argv[]) {
	//args = input.fastq output_dir spacer.fasta fragment_size tag_size
	if(argc!=6){
		printf("usage : reconstruction.cpp input.fastq output_dir spacer.fasta fragment_size tag_size");
		return 1;
	}


	string input_path = argv[1];
	string output_dir = argv[2];
	string spacer_path = argv[3];

	string SPACER = read_spacer(spacer_path);
	vector<string> SEQ = read_sequences(input_path);
	int FRAG_SIZE = stoi(argv[4]);
	int TAG_SIZE = stoi(argv[5]);

	printf("Spacer:  %s\n", SPACER.c_str());
	printf("Fragment size:  %d\n", FRAG_SIZE);
	printf("Tag size:  %d\n", TAG_SIZE);

	reconstruction(SEQ, SPACER, FRAG_SIZE, output_dir);

	printf("done !");

	return 0;
}

string read_spacer(string spacer_path) {
	ifstream file(spacer_path);
	string line;
	getline(file, line);
	getline(file, line);
	file.close();
	return line;
}

vector<string> read_sequences(string sequences_path) {
	//read the sequences from a fastq file
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

string* add_kmer_to_dict(string &kmer) {
	if(kmer_dict[kmer] == NULL) {
		string* p =(string*) malloc(sizeof(kmer));
		kmer_dict[kmer] = p;
		//cout << "added kmer " << *kmer << "->"<< p << endl;
		return p;
	} else {
		//cout << "kmer already in dict " << *kmer << "->"<< kmer_dict[*kmer] << endl;
		return kmer_dict[kmer];
	}
}

void reconstruction(vector<string> SEQ, string SPACER, int FRAG_size, string output_dir) {

	//CUT SEQUENCES
	auto start = chrono::system_clock::now();


	int WS = int(SPACER.size())/2; //size of the kmers

	set<string> K_spacer;
	set<string> K_rev_spacer;
	string RC_SPACER = reverse_complement(SPACER);
	for(int k=0; k<WS+1; k++) {
		string kmer = SPACER.substr(k, WS);
		K_spacer.insert(kmer);
		kmer = RC_SPACER.substr(k, WS);
		K_rev_spacer.insert(kmer);
	}
	vector<string> FRAG;
	for(const string& sequence : SEQ) {
		// locate spacer words
		vector<tuple<int, string>> L;
		int end_range = sequence.size()-WS+1;
		for(int j=0; j<end_range; j++) {
			string kmer = sequence.substr(j, WS);
			if(K_spacer.find(kmer) != K_spacer.end()
					or K_rev_spacer.find(kmer) != K_rev_spacer.end()) {
				L.push_back(make_tuple(j, kmer));
			}
		}
		// group close spacer coordinates together
		vector<vector<tuple<int, string>>> COORD;
		int last_coord = -100;
		vector<tuple<int, string>> set_coord;
		for(const tuple<int, string>& couple : L) {
			int coord = get<0>(couple);
			string kmer = get<1>(couple);
			if(coord - last_coord < int(SPACER.size()/2)+1) {
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
			FRAG.push_back(SPACER+fragment+SPACER); //re add the original spacer for the consensus
		}
	}
	auto end = chrono::system_clock::now();

	chrono::duration<double> elapsed_seconds = end-start;
	time_t end_time = chrono::system_clock::to_time_t(end);

	cout << "cut sequences: " << elapsed_seconds.count() << "s\n";
	//ASSIGN SET KMER
	start = chrono::system_clock::now();


	cout << "assign_set_kmer " << FRAG.size() << endl;

	int len_spacer = SPACER.size();
	map<string*, int> kmers_p;
	for(const string& fragment : FRAG) {
		int end_range = fragment.size()-len_spacer-WS+1;
		for(int j=len_spacer; j<end_range; j++) {
			string kmer = fragment.substr(j, WS);
			string* p = add_kmer_to_dict(kmer);
			kmers_p[p]++;
		}
	}
	int SOLID_THRESHOLD = 5;
	//assign a set of kmer to each fragment
	vector<vector<string*>> list_kmer_p;
	for(const string& fragment : FRAG) {
		set<string*> frag_set;
		int end_range = fragment.size()-len_spacer-WS+1;
		for(int j=len_spacer; j<end_range; j++) {
			string kmer = fragment.substr(j, WS);
			string* p = kmer_dict[kmer];
			if(kmers_p[p] >= SOLID_THRESHOLD){
				frag_set.insert(p);
			}
		}
		vector<string*> frag_vector(frag_set.size());
		copy(frag_set.begin(), frag_set.end(), frag_vector.begin());
		list_kmer_p.push_back(frag_vector);
	}

	end = chrono::system_clock::now();

	elapsed_seconds = end-start;
	end_time = chrono::system_clock::to_time_t(end);

	cout << "assign kmers: " << elapsed_seconds.count() << "s\n";


	//CONSTRUCT FRAGMENT GRAPH

	start = chrono::system_clock::now();

	int len_FRAG = FRAG.size();
	cout << "construct_fragment_graph " << list_kmer_p.size() << endl;
	Graph G(len_FRAG);
	//int progress = -1;
	for(int i=0; i<len_FRAG; i++) {
		/*if(i*100/len_FRAG > progress) {
			progress = i*100/len_FRAG;
			printf(" %d %%\n",progress);
		}*/
		int size_a = list_kmer_p[i].size();

		for(int j=i+1; j<len_FRAG; j++) {
			//int intersection_size = count_shared_elements(set_kmer_p[i], set_kmer_p[j]);
			int size_b = list_kmer_p[j].size();


			//return the number of shared elements between 2 sorted vectors of pointers
			int index_a = 0;
			int index_b = 0;
			int intersection_size = 0;
			while(index_a < size_a and index_b < size_b) {
				//cout << vector_a[index_a] << " " << vector_b[index_b] << endl;
				string* l_a = list_kmer_p[i][index_a];
				string* l_b = list_kmer_p[j][index_b];
				if(l_a < l_b) {
					index_a++;
				} else if(l_a > l_b) {
					index_b++;
				} else {
					index_a++;
					index_b++;
					intersection_size++;
					//cout << vector_a[index_a] << " " << vector_b[index_b] << endl;
				}
			}

			//add an edge in the graph if the 2 fragments share enough kmers
			if(intersection_size > FRAG_size/3){
				//cout <<"size : "<< intersection_size << endl;
				G.add_edge(i, j);
			}
		}
	}

	end = chrono::system_clock::now();

	elapsed_seconds = end-start;
	end_time = chrono::system_clock::to_time_t(end);

	cout << "construct frag graph: " << elapsed_seconds.count() << "s\n";

	//PRINT CONNECTED COMPONENTS

	//save the list of clusters
	cout << "print_connected_components" << endl;
	int k = 0;
    cout << "mkdir "<< output_dir+"/frag : " << mkdir((output_dir+"/frag").c_str(), 0777) << endl;
    vector<vector<int>> connected_components = G.connected_components(G.v);
	cout << "connected components " << connected_components.size() << endl;
    for(const vector<int>& cc : connected_components) {
    	//cout << cc.size() << endl;
    	if(cc.size() > 10) {
    		ofstream outfile(output_dir+"/frag/f"+to_string(k)+".fasta");
    		for(const int& node : cc) {
    			//cout << node << ".";
    			outfile << ">s"+to_string(k)+"_"+to_string(node)+"\n"+FRAG[node] << endl;
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
