#include "weisfeiler_lehman.h"
#include <vector>
#include <iostream>
#include <math>
namespace wl{
        private class AdjMatrix{
                public:
                        AdjMatrix(int size):num_of_vertices(size){
                                int numberOfBits = size * size;
                                int numOfWords = std::ceil(numberOfBits / (8*sizeof(int)));
                                array = new int[numOfWords];
                        }
                        ~AdjMatrix(){
                                del[] array;
                        }
                        void addEdge(int v, int u, bool bothDirections=false){
                                unsigned int bit = getEdgeBit(v,u);
                                int word, wordOffset;
                                wordOffset = bit % wordBits;
                                word = bit / wordBits;
                                unsigned int mask = (1<<wordOffset);
                                array[word] |= mask;
                                if(bothDirections) addEdge(u,v);
                        }
                        void delEdge(int v, int u, bool bothDirections=false){
                                unsigned int bit = getEdgeBit(v,u);
                                int word, wordOffset;
                                wordOffset = bit % wordBits;
                                word = bit / wordBits;
                                unsigned int mask = -1;
                                mask ^= (1<<wordOffset);
                                array[word] &= mask;
                                if(bothDirections) delEdge(u,v);
                        }
                        bool isEdge(int v, int u){
                                unsigned int bit = getEdgeBit(v,u);
                                int word, wordOffset;
                                wordOffset = bit % wordBits;
                                word = bit / wordBits;
                                unsigned int mask = (1<<wordOffset);
                                unsigned int res = array[word] & mask;
                                return res != 0;
                        }
                private:
                        int num_of_vertices;
                        int* array;
                        inline unsigned int getEdgeBit(int v, int u){
                                int v_row_offset = v * num_of_vertices;
                                return v_row_offset + u;
                        }
                        const int wordBits = (8*sizeof(int));
        };
        
        int prova2(std::vector<GraphNode> v){
                for(const auto& el: v){
                        std::cout << el.idx << " " << el.color << "->";
                        for(const auto& el2: el.adj_list) std::cout << el2 << " ";                        
                        std::cout << std::endl;
                }
                return v.size();
        }
}
