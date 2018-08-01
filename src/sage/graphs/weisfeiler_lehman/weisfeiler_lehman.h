#ifndef WL_H
#define WL_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include "Tuple.h"
namespace wl{
        struct GraphNode{
                long long idx, color;
                std::vector<pair<int,int>> adj_list;
        };
        
        struct IntPair_Hash {
            size_t operator() (const pair<int,int>& p) const {
                size_t h = 0;
                hash_combine(h, p.first);
                hash_combine(h, p.second);
                return h;
            }
        };

        unordered_map<int, vector<pair<int,int>>> k_WL(const std::vector<GraphNode>& v, int k, bool hasVertexLabels = false);
        struct IntVector_Hash {
            size_t operator() (const vector<int>& p) const {
                size_t h = 0;
                for(const auto& el: p){
                    hash_combine(h, el);
                }
                return h;
            }
        };
}
#endif
