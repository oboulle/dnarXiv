#include <iostream>
#include <list>
using namespace std;


class Graph{

	public:

    int size;
    list<int> *edges;

    Graph(int V){
        size=V;
        edges=new list<int> [size];
    }

    void add_edge(int u,int v){
        edges[u].push_back(v);
        edges[v].push_back(u);
    }

    void print(){
        for(int node=0; node<size; node++){
            cout<< node <<" -> ";
            for(auto neighbour=edges[node].begin(); neighbour!=edges[node].end(); neighbour++){
                cout<< *neighbour <<" ";
            }
            cout<<"\n";
        }
    }

    void get_neighbourhood(int node, bool visited[], vector<int>& neighbourhood){
    	neighbourhood.push_back(node);
        visited[node] = true;
        for(auto neighbour=edges[node].begin(); neighbour!=edges[node].end(); neighbour++){
            if(!visited[*neighbour]){
            	get_neighbourhood(*neighbour, visited, neighbourhood);
            }
        }
    }

    vector<vector<int>> get_connected_components(){
    	//get a list of lists of connected components in the graph
       	vector<vector<int>> connected_components;
       	bool* visited = new bool[size]();
        for(int node=0; node<size; node++){
            if(!visited[node]){
            	vector<int> neighbourhood;
            	get_neighbourhood(node, visited, neighbourhood);
            	connected_components.push_back(neighbourhood);
            }
        }
        return connected_components;
    }

    vector<vector<int>> filter_low_connected_components(vector<vector<int>> connected_components){
    	//remove the nodes that have a weak connection from a list of connected components
    	vector<vector<int>> filtered_components;
    	float edges_mean = get_edges_mean();
    	for(int i=0; i<connected_components.size(); i++){
    		vector<int> filtered_nodes;
    		for(const int& node : connected_components[i]) {
    			//weak connection = low number of edges compared to the mean
    			if (edges[node].size() > edges_mean/2.0) {
    				filtered_nodes.push_back(node);
    			}
    		}
    		//filter empty components
    		if (filtered_nodes.size() > 0) {
    			filtered_components.push_back(filtered_nodes);
    		}
    	}
    	return filtered_components;
    }

    float get_edges_mean() {
    	int sum=0;
    	for(int node=0; node<size; node++){
    		sum = sum + edges[node].size();
    	}
    	return float(sum)/float(size);
    }
};


int main_test()
{
    int v=6;
    Graph g(v);
    g.add_edge(0,1);
    g.add_edge(1,2);
    g.add_edge(2,3);
    g.add_edge(3,1);
    g.add_edge(4,5);

    return 0;
}
