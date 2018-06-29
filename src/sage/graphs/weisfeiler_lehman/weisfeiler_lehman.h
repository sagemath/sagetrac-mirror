#ifndef WL_H
#define WL_H

#include<vector>

namespace wl{
        struct GraphNode{
                long long idx, color;
                std::vector<int> adj_list;
        };
        bool k_WL(const std::vector<GraphNode>& v, int k);
}
#endif
