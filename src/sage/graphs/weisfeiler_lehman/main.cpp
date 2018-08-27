#include <fstream>
#include <vector>
#include "weisfeiler_lehman.h"
#include <unordered_set>
using namespace std;
int main(int argc, char** argv){
	ifstream graph("input");
	int k, card;
	if(argc == 1){
		cout << "Insert k value: ";
		cin >> k;
        cout << "Insert cardinality: ";
        cin >> card;
	}else{
		k = atoi(argv[1]);
        card = atoi(argv[2]);
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
	auto res = wl::k_WL(v, k, false);
    auto res2 = res;
    if(card == 1){
        res2.clear();
        for(const auto& el: res){
            vector<pair<int,int>> temp;
            for(const auto& el2: el.second){
                if(el2.first == el2.second) temp.push_back(el2);
            }
            if(!temp.empty()) res2[el.first] = temp;
        }
    }
    /*for(const auto& el: res2){
        cout << el.first << ": " << endl;
        for(const auto& el2: el.second){
            cout << "       (" << el2.first << ", " << el2.second << ")" << endl;
        }
    }*/
    cout << res2.size() << " classi di colore" << endl;
	return 0;
}
