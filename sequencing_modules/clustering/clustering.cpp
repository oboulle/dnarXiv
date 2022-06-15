//============================================================================
// Name        : Clustering.cpp
// Author      : Olivier
//============================================================================

//g++ -std=c++11 clustering.cpp -o clustering


#include <map>
#include <algorithm>
#include <vector>
#include "graph.h"
#include <fstream>
#include <set>
#include <chrono>

using namespace std;

string read_spacer(string spacer_path);
vector<string> read_sequences(string sequences_path);
vector<string> fragments_identification(vector<string> sequences, string spacer, int fragment_size);
void clustering(vector<string> fragments_vector, string spacer, int fragment_size, string output_frag_dir);


int main(int argc,char* argv[]) {
	//args = input.fastq output_frag_dir spacer.fasta fragment_size
	if(argc!=5){
		printf("usage : clustering input.fastq output_frag_dir spacer fragment_size\n");
		return 1;
	}
	printf("clustering...\n");
	string input_path = argv[1]; //fastq file of the input sequences
	string output_frag_dir = argv[2]; //directory to save the cluster files
	string spacer = argv[3]; //spacer sequence
	int fragment_size = stoi(argv[4]); //original size of the fragments (without spacers)
	vector<string> sequences = read_sequences(input_path);

	//printf("Spacer:  %s\n", spacer.c_str());
	//printf("Fragment size:  %d\n", fragment_size);

	//vector<string> sequences = vector<string>(sequences.begin(), sequences.begin() + sequences.size()/2); // @suppress("Symbol is not resolved")
	//vector<string> sequences_2 = vector<string>(sequences.begin() + sequences.size()/2, sequences.end()); // @suppress("Symbol is not resolved")
	//REMOVED//apply a clustering to the first half of the sequences and another independent clustering to the second half

	vector<string> fragments_vector = fragments_identification(sequences, spacer, fragment_size);
	clustering(fragments_vector, spacer, fragment_size, output_frag_dir);
	//clustering(fragments_vector_2, spacer, fragment_size, output_frag_dir");
	printf("\tcompleted !\n");
	return 0;
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
	for (size_t i = 0; i < sequence.length(); ++i) {
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


vector<string> fragments_identification(vector<string> sequences, string spacer, int fragment_size) {
	//cut the fragments from the sequences
	//return a list of fragments
	auto start = chrono::system_clock::now();

	int len_spacer = spacer.size();
	int len_kmer = int(len_spacer/2); //length of the kmers

	string rev_spacer = reverse_complement(spacer);

	set<string> spacer_kmers; //set of kmers from the spacer
	set<string> rev_spacer_kmers; //set of kmers from the reverse complement spacer

	//init the sets of kmers from the spacer and the reverse spacer
	for(int i=0; i<len_spacer-len_kmer+1; i++) {
		string kmer = spacer.substr(i, len_kmer);
		spacer_kmers.insert(kmer);
		kmer = rev_spacer.substr(i, len_kmer);
		rev_spacer_kmers.insert(kmer);
	}
	vector<string> fragments_vector;
	for(const string& sequence : sequences) {
		// locate coordinates of (rev)spacer kmers
		vector<tuple<int, string>> spacer_kmers_coordinates; //contains pairs (coord;kmer) of spacers kmers in the sequence
		int end_range = sequence.size()-len_kmer+1;
		for(int i=0; i<end_range; i++) {
			string kmer = sequence.substr(i, len_kmer);
			if(spacer_kmers.find(kmer) != spacer_kmers.end()
					or rev_spacer_kmers.find(kmer) != rev_spacer_kmers.end()) {
				spacer_kmers_coordinates.push_back(make_tuple(i, kmer));
			}
		}
		// group close (rev)spacer kmers coordinates together
		vector<tuple<int, string>> close_couples; //vector of couples (coord/kmer) that are close to each other
		vector<vector<tuple<int, string>>> COORD; //vector of close couples vectors

		int last_coord = -1000;
		for(const tuple<int, string>& couple_coord_kmer : spacer_kmers_coordinates) {
			int coord = get<0>(couple_coord_kmer);
			string kmer = get<1>(couple_coord_kmer);
			if(coord - last_coord < len_kmer+1) { //the couple is located close to the previous couple
				close_couples.push_back(couple_coord_kmer); //add in the current vector of close couples
			} else {
				if(close_couples.size() != 0) {
					//add the previous vector of close couples to COORD
					int l = COORD.size();
					COORD.push_back(vector<tuple<int, string>>());
					for(const tuple<int, string>& couple : close_couples) {
						COORD[l].push_back(couple);
					}
				}
				close_couples = vector<tuple<int, string>> {couple_coord_kmer}; //init a new vector of close couples
			}
			last_coord = coord;
		}
		//add the last vector of close couples to COORD
		int l = COORD.size();
		COORD.push_back(vector<tuple<int, string>>());
		for(const tuple<int, string>& couple : close_couples) {
			COORD[l].push_back(couple);
		}

		//determine spacer location
		for(int i=0; i<COORD.size()-1; i++) {

			tuple<int, string> end_spacer = COORD[i].back(); //spacer before the fragment
			tuple<int, string> start_next_spacer = COORD[i+1][0]; //spacer after the fragment


			int end_spacer_coord = get<0>(end_spacer) + get<1>(end_spacer).size(); //coord where the spacer before the fragment stop
			int start_next_spacer_coord = get<0>(start_next_spacer); //coord where the spacer after the fragment starts

			string fragment = sequence.substr(end_spacer_coord,start_next_spacer_coord-end_spacer_coord);
			//cout << fragment << "\n" ;

			if(fragment.size() > int(2*fragment_size/3) and fragment.size() < int(3*fragment_size/2)) { //ignore fragments that are too small (non existant spacer recognized) or too large (spacer not recognized between 2 fragments)

				if(rev_spacer_kmers.find(get<1>(end_spacer)) != rev_spacer_kmers.end()
										and rev_spacer_kmers.find(get<1>(start_next_spacer)) != rev_spacer_kmers.end()) {
					//the fragment has been reverted
					fragment = reverse_complement(fragment);
				}
				fragments_vector.push_back(fragment);
			}
		}
	}
	chrono::duration<double> elapsed_seconds = chrono::system_clock::now()-start;
	//cout << "cut sequences: " << elapsed_seconds.count() << "s\n";

	return fragments_vector;
}


void clustering(vector<string> fragments_vector, string spacer, int fragment_size, string output_frag_dir) {
	//ASSIGN SET KMER
	//assign a set of kmer to each fragment
	auto start = chrono::system_clock::now();
	//cout << "assign_set_kmer " << fragments_vector.size() << endl;

	int len_kmer = int(spacer.size()/2); //length of the kmers

	map<string, int> list_kmers; //map of pairs kmers;nbr_occurrence in all the fragments

	for(const string& fragment : fragments_vector) {
		int end_range = fragment.size()-len_kmer+1;
		string kmer;
		for(int i=0; i<end_range; i++) {
			kmer = fragment.substr(i, len_kmer);
			list_kmers[kmer]++;
		}
	}
	int min_kmers_threshold = 5; //minimum number of occurrence for a kmer in list_kmers

	vector<vector<string>> set_kmers_per_fragment; //list of kmers set, index = fragment index in the list of fragments
	for(const string& fragment : fragments_vector) {
		vector<string> related_kmers; //kmers of the fragment
		int end_range = fragment.size()-len_kmer+1;
		string kmer;
		for(int i=0; i<end_range; i++) {
			kmer = fragment.substr(i, len_kmer);
			if(list_kmers[kmer] >= min_kmers_threshold) {
				related_kmers.push_back(kmer);
			}
		}
		set_kmers_per_fragment.push_back(related_kmers);
	}

	chrono::duration<double> elapsed_seconds = chrono::system_clock::now()-start;
	//cout << "assign kmers: " << elapsed_seconds.count() << "s\n";

	//CONSTRUCT FRAGMENT GRAPH / clustering
	start = chrono::system_clock::now();

	int len_frag = fragments_vector.size();
	//cout << "construct_fragment_graph " << set_kmers_per_fragment.size() << endl;
	Graph G(len_frag); //init a graph, node i = fragment of index i
	map<string, vector<int>> kmer_occurence; //map of pairs kmer;(list of index for fragments containing this kmer)

	for(int i=0; i<len_frag; i++) {
		for(const string& kmer : set_kmers_per_fragment[i]) {
			kmer_occurence[kmer].push_back(i);
		}
	}

	//threshold for number of shared kmers between 2 fragments to be linked by an edge
	//a fragment has =~ <fragment_size> overlapping kmers
	float shared_kmers_rate = (1.0/3.0);
	for(int i=0; i<len_frag; i++) {
		vector<int> shared_kmers(len_frag, 0); //vector where shared_kmers[i] = number of shared kmers with fragment i
		for(const string& kmer : set_kmers_per_fragment[i]) {
			for(const int& num_frag : kmer_occurence[kmer]) {
				shared_kmers[num_frag]++;
			}
		}
		for(int num_frag=i+1; num_frag<len_frag; num_frag++) {
			if(shared_kmers[num_frag] > fragment_size*shared_kmers_rate) {
				//an edge is created between 2 fragments if enough kmers are shared
				G.add_edge(i, num_frag);
			}
		}
	}

	elapsed_seconds = chrono::system_clock::now()-start;
	//cout << "construct frag graph: " << elapsed_seconds.count() << "s\n";
	start = chrono::system_clock::now();

	//PRINT CONNECTED COMPONENTS
	//save the list of clusters fragments
	int k = 0;
	//cout << "mkdir "<< output_frag_dir+" : " << mkdir((output_frag_dir).c_str(), 0777) << endl;
    vector<vector<int>> connected_components = G.get_connected_components();

    vector<vector<int>> filtered_components = connected_components;//G.filter_low_connected_components(connected_components);
	//cout << "connected components " << connected_components.size() << endl;
    for(const vector<int>& cluster : filtered_components) {
    	//cout << cc.size() << endl;
    	if(cluster.size() > 10) {
    		ofstream outfile(output_frag_dir+"/f"+to_string(k)+".fasta");
    		for(const int& node : cluster) {
    			outfile << ">s"+to_string(k)+"_"+to_string(node)+"\n"+fragments_vector[node] << endl;
    		}
    		//cout << endl;
    		outfile.close();
    		k++;
    	}
    }
	//cout << k << " clusters" << endl;
    elapsed_seconds = chrono::system_clock::now()-start;
    //cout << "save connected components: " << elapsed_seconds.count() << "s\n";
}
