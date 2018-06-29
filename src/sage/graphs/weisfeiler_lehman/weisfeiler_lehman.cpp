#include "weisfeiler_lehman.h"
#include "Tuple.h"
#include <vector>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <cstring>
namespace wl{
        using TupleMap = std::unordered_map<Tuple<int>, int>;
        using InverseTupleMap = std::vector<std::unordered_map<Tuple<int>, int>::const_iterator>;
        using uchar = unsigned char;
        class AdjMatrix{
                public:
                        AdjMatrix(int size):num_of_vertices(size){
                                int numberOfBits = size * size;
                                int numOfWords = std::ceil(numberOfBits / (8*sizeof(int)));
                                array = new int[numOfWords];
                        }
                        ~AdjMatrix(){
                                delete[] array;
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
                        bool isEdge(int v, int u) const{
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
                        inline unsigned int getEdgeBit(int v, int u) const{
                                int v_row_offset = v * num_of_vertices;
                                return v_row_offset + u;
                        }
                        const int wordBits = (8*sizeof(int));
        };
        class AtomicType{
                public:
                        template<typename Container>
                        AtomicType(const Container& vertex_subset, const AdjMatrix& adjMatrix, int originalColor = 0){
                              this->vertex_subset.reserve(vertex_subset.size());
                              std::copy(vertex_subset.begin(), vertex_subset.end(), this->vertex_subset.begin());
                              this->size = vertex_subset.size();
                              this->original_color = originalColor;
                              atp_matrix.resize(this->size*this->size);
                              for(int i = 0; i < this->size; i++){
                                      for( int j = 0; j < this->size; j++){
                                                int iV = this->vertex_subset[i];
                                                int jV = this->vertex_subset[j];
                                                int atpIdx = getIdx(iV,jV);
                                                if(iV == jV) atp_matrix[atpIdx] = 2;
                                                else atp_matrix[atpIdx] = adjMatrix.isEdge(iV, jV)?1:0;
                                      }
                              }
                        }
                        bool operator==(const AtomicType& b) const{
                                if(size != b.size) return false;
                                if(original_color != b.original_color) return false;
                                return atp_matrix == b.atp_matrix;
                        }
                        const vector<uchar>& getSetVector() const{
                                return atp_matrix;
                        }
                private:
                        inline int getIdx(int row, int column){
                                return row*size + column;
                        }
                        int size;
                        int original_color;
                        vector<int> vertex_subset;
                        vector<uchar> atp_matrix;
        };
        class Coloring{
                public:
                        Coloring(int pc): previous_color(pc), set_hash(0){};
                        void update_previous_color(int newColor){
                                set_coloring.clear();
                                set_hash = 0;
                                previous_color = newColor;
                        }
                        void add_multiset_color(int color){
                                set_coloring.push_back(color);
                                hash_combine<int>(set_hash, color, true);
                        }
                        bool operator==(const Coloring& b) const{
                                if(previous_color != b.previous_color) return false;
                                if(set_hash != b.set_hash) return false;
                                return set_coloring == b.set_coloring;
                        }
                        const vector<int>& getSetVector() const{
                                return set_coloring;
                        }
                private:
                        int previous_color;
                        size_t set_hash;
                        vector<int> set_coloring;
        };
        template<typename T>
        unordered_map<int, unordered_set<int>> createEquivalenceClasses(const vector<T>& whole_set){
                unordered_map<int, unordered_set<int>> result;
                unordered_set<int> remaining;
                int size = whole_set.size();
                for(int i = 0; i < size; i++){
                        remaining.insert(i);
                }
                while(!remaining.empty()){
                        int representative = *(remaining.begin());
                        remaining.erase(representative);
                        const T& repr_element = whole_set.at(representative);
                        result[representative];
                        for(const auto& idx: remaining){
                                const T& set_el = whole_set.at(idx);
                                if(set_el == repr_element) result[representative].insert(idx);
                        }
                        for(const auto& idx: result[representative]){
                                remaining.erase(idx);
                        }
                }
                return result;
        }

        void innerGenerateTupleMap(int i, int k, int n, vector<int>& currentTuple, TupleMap& tupleMap, InverseTupleMap& inverseTupleMap){
                if(i == k){
                        auto t = Tuple<int>(currentTuple);
                        size_t idx = tupleMap.size();
                        auto tuple_iterator = tupleMap.emplace(t, idx).first;
                        inverseTupleMap[idx] = tuple_iterator;
                        return;
                }
                for(int j = 0; j < n; j++){
                        currentTuple[i] = j;
                        innerGenerateTupleMap(i+1, k, n, currentTuple, tupleMap, inverseTupleMap);
                }
                return;
        }
        
        void generateTupleMap(int n, int k, TupleMap& tupleMap, InverseTupleMap& inverseTupleMap){
                inverseTupleMap.resize((size_t)(pow(n,k)));
                vector<int> tempTuple;
                tempTuple.resize(k);
                innerGenerateTupleMap(0, k, n, tempTuple, tupleMap, inverseTupleMap);
        }
        
        class prova{
                vector<char> v;
                public:
                const vector<char>& getSetVector() const{
                        return v;
                }
                prova(int s){
                        v.resize(s);
                        for(auto& el: v){
                                el = rand()% 2;
                        }
                }
                prova(){
                        v.resize(10000);
                        for(auto& el: v){
                                el = rand()% 100;
                        }
                }
                bool operator<(const prova& b) const{
                        for(size_t i = 0; i < v.size(); i++){
                                if(v[i] == b.v[i]) continue;
                                return v[i] < b.v[i];
                        }
                        return false;
                }
                bool operator==(const prova& b) const{
                        return v == b.v;
                }
        };
        std::ostream &operator<<(std::ostream &stream, const prova& t) {
                        for(const auto& el:t.getSetVector()){
                                stream << el << ", ";
                        }
                        stream << std::endl;
                        return stream;
        }
        
        template<typename T>
        vector<vector<int>> countingSort(vector<pair<T,int>>& toOrder, T m, T M){
                std::cerr<<"counting" << std::endl;
                vector<vector<int>> buckets(M-m+1);
                int s = toOrder.size();
                for(int i = 0; i < s; i++){
                        auto& element = toOrder[i];
                        buckets[element.first-m].push_back(element.second);
                }
                return buckets;
        }
        template<typename T>
        vector<vector<int>> generateBuckets(vector<pair<T,int>>& valuesToBucket){
                //std::cerr<<"quick" << std::endl;
                vector<vector<int>> res;
                T lastValue = 0;
                for(auto& el: valuesToBucket){
                        if(res.size() == 0 || lastValue != el.first){
                                res.push_back(std::vector<int>{});
                        }
                        res.back().push_back(el.second);
                        lastValue = el.first;
                }
                return res;
        } 
        template<class T>
        vector<vector<int>> orderSliceByElement(const vector<T>& sortedSets, vector<int>::const_iterator mapSliceB, vector<int>::const_iterator mapSliceE, int sortingIndex){
                auto minElement = sortedSets[*mapSliceB].getSetVector()[sortingIndex], maxElement = sortedSets[*mapSliceB].getSetVector()[sortingIndex];
                vector<pair<decltype(minElement),int>> elementsToOrder; //First value is the value of the orderingIndex-th element, second value is the index of the set in sortedSets
                //This copy could be removed as an optimization, but at the moment it's not worth the effort. In case, parameters for minElement and maxElement will be needed
                for(auto it = mapSliceB; it != mapSliceE; it++){
                        auto element = sortedSets[*it].getSetVector()[sortingIndex];
                        if(minElement > element) minElement = element;
                        if(maxElement < element) maxElement = element;
                        elementsToOrder.push_back({element, *it});
                }
                auto s = elementsToOrder.size();
                vector<vector<int>> res;
                if(s * log2(s) < maxElement-minElement || !(std::is_convertible<T,int>::value)){
                        sort(elementsToOrder.begin(), elementsToOrder.end());
                        res = generateBuckets(elementsToOrder);
                }else{
                        res = countingSort(elementsToOrder, minElement, maxElement);
                }
                return res;
        }
        template<typename T>
        void iota(typename std::vector<T>::iterator b, typename std::vector<T>::iterator e, int start){
                for(auto it = b; it != e; it++){
                        *it = start++;
                }
        }
        template<class T>
        vector<int> orderSortedSets(const vector<T>& sortedSets, int setSize){
                using BucketTuple = std::tuple<int,int,int>;
                vector<int> orderMap(sortedSets.size());
                iota<int>(orderMap.begin(), orderMap.end(), 0);
                queue<BucketTuple> bucketsToSort;
                bucketsToSort.push({0, sortedSets.size(), 0});
                while(!bucketsToSort.empty()){
                        BucketTuple v = bucketsToSort.front();
                        int tupleB = std::get<0>(v);
                        int tupleE = std::get<1>(v);
                        int sortingIndex = std::get<2>(v);
                        //std::cout << "Bucket: (" << tupleB << ", " << tupleE << ", " << sortingIndex << ")" << std::endl;
                        bucketsToSort.pop();
                        auto res = orderSliceByElement(sortedSets, orderMap.cbegin()+tupleB, orderMap.cbegin()+tupleE, sortingIndex);
                        int b = 0;
                        for(auto& bucket: res){
                                int bucketSize = bucket.size();
                                if(bucketSize == 0) continue;
                                int e = b+bucketSize;
                                std::swap_ranges(orderMap.begin()+tupleB+b, orderMap.begin()+tupleB+e, bucket.begin());
                                //std::cout << "------- SottoBucket: (" << tupleB+b << ", " << tupleB+e << ", " << sortingIndex+1 << ")" << std::endl;
                                if(e - b > 1 && sortingIndex + 1 != setSize) bucketsToSort.push({tupleB+b, tupleB+e, sortingIndex+1});
                                b = e;
                        }
                }
                return orderMap;
        }
        
        bool k_WL(const std::vector<GraphNode>& v, int k){
                AdjMatrix adj_matrix(v.size());
                for(const auto& el: v){
                        int vIdx = el.idx;
                        for(const auto& adj: el.adj_list){
                                adj_matrix.addEdge(vIdx, adj);
                        }
                }
                TupleMap tm;
                InverseTupleMap itm;
                generateTupleMap(v.size(), k, tm, itm);
                vector<AtomicType> atp;
                atp.reserve(itm.size());
                for(const auto& el: itm){
                        atp.push_back(AtomicType(el->first, adj_matrix));
                }
                vector<int> atp_remap = orderSortedSets(atp, k*k);
                
                //Now the idea is going through the atps in remapping order,
                //and initialize the color of each tuple based on the index of the bucket the tuple is in
                
                //After this, the main part of the algorithm, comprised of computing the coloring of adjacents, ordering the
                //tuples by old_coloring^multiset and then update their old coloring.
                //One should stop when the orbits are the same before and after a round. This last part is gonna be tricky for sure
                
                
                /*for(const auto& el: itm){
                        std::cout << el->second << " = " << el->first << std::endl;
                }*/
                vector<prova> v2(100);
                vector<prova> v3 = v2;
                /*for(const auto& el: v2){
                        std::cout << el << std::endl;
                }*/
                sort(v2.begin(), v2.end());
                vector<int> remap = orderSortedSets(v3, 10000);
                for(size_t i = 0; i < v2.size(); i++){
                        if(v2[i].getSetVector() != v3[remap[i]].getSetVector()) return false;
                        /*        std::cout << "Il numero " << i << " non è uguale" << std::endl;
                                std::cout << " --------- " << v2[i] << std::endl;
                                std::cout << " ***************************************************************************************";
                                std::cout << " --------- " << v3[remap[i]] << std::endl;
                        }else std::cout << "Il numero " << i << " è uguale" << std::endl;*/
                }
                return true;
        }

}
