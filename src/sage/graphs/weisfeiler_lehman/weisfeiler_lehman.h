#ifndef WL_H
#define WL_H

#include <vector>
#include <unordered_map>
#include "Tuple.h"
namespace wl{
        using TupleMap = std::unordered_map<Tuple<int>, int>;
        struct GraphNode{
                long long idx, color;
                std::vector<int> adj_list;
        };
        TupleMap k_WL(const std::vector<GraphNode>& v, int k);
}
#endif
