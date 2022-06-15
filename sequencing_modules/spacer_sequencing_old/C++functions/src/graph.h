#include <iostream>
#include <list>
using namespace std;

using namespace std;

class Graph{
public:
    int v;
    list<int> *l;
    Graph(int V){
        v=V;
        l=new list<int> [v];
    }
    void add_edge(int u,int v){
        l[u].push_back(v);
        l[v].push_back(u);
    }
    void print(){
        for(int i=0;i<v;i++){
            cout<<i<<" -> ";
            for(auto it=l[i].begin();it!=l[i].end();it++){
                cout<<*it<<" ";
            }
            cout<<"\n";
        }
    }
    void dfsUtil(int src,bool visited[],vector<int>& s){
        s.push_back(src);
        visited[src]=true;
        for(auto it=l[src].begin();it!=l[src].end();it++){
            if(!visited[*it]){
                dfsUtil(*it,visited,s);
            }
        }
    }

    vector<vector<int>> connected_components(int v){
       	vector<vector<int>> cc;
       	bool* visited = new bool[v]();
        //int c=1;
        for(int i=0;i<v;i++){
            if(!visited[i]){
            	vector<int>s;
            	dfsUtil(i,visited,s);
            	cc.push_back(s);
                //cout<<"Component "<< c<<" : " ;
                //for(auto j=s.begin();j!=s.end();j++){
                //    cout<<*j<<",";
                //}
                //cout<<"\n";
                //c++;
           }
        }
        return cc;
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
    g.connected_components(v);
    return 0;
}
