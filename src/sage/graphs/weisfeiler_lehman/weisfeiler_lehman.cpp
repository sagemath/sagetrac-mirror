#include "weisfeiler_lehman.h"
#include <vector>
#include <iostream>
namespace wl{
        int prova2(std::vector<GraphNode> v){
                for(const auto& el: v){
                        std::cout << el.idx << " " << el.color << "->";
                        for(const auto& el2: el.adj_list) std::cout << el2 << " ";                        
                        std::cout << std::endl;
                }
                return v.size();
        }
}
