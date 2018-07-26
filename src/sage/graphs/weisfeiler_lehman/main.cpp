#include <fstream>
#include <vector>
#include "weisfeiler_lehman.h"
using namespace std;
int main(int argc, char** argv){
	ifstream graph("input");
	int k;
	if(argc == 1){
		cout << "Insert k value: ";
		cin >> k;
	}else{
		k = atoi(argv[1]);
	}
	vector<wl::GraphNode> v;
	int n, colors;
	graph >> colors >> n;
	for(int i = 0; i < n; i++){
		wl::GraphNode g;
		g.idx = i;
		g.color = 1;
		for(int j = 0;j < n; j++){
			int edge;
			graph >> edge;
			if(i == j) g.color = edge+1;
			else if(edge != 0) g.adj_list.push_back(make_pair(j,edge));
		}
		v.push_back(g);
	}
	auto res = wl::k_WL(v, k);
	#ifdef DEBUG
	for(const auto& el: res){
		cout << el.first << ":" << endl;
        for(const auto& el2: el.second){
            cout << "    (" << el2.first << ", " << el2.second << ")" << endl;
        }
	}
	#endif
	return 0;
}
