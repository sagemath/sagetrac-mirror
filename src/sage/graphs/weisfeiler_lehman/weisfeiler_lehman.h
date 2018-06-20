#ifndef WL_H
#define WL_H

#include<vector>

namespace wl{
        struct GraphNode{
                long long idx, color;
                std::vector<int> adj_list;
        };
        int prova2(std::vector<GraphNode> v);
}
#endif
