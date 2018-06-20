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
#include <functional>
namespace wl{
        using TupleMap = std::unordered_map<Tuple<int>, int>;
        using InverseTupleMap = std::vector<std::unordered_map<Tuple<int>, int>::const_iterator>;
        
        template<typename T>
        class AdjMatrix{
                public:
                        //The type T default constructor should set an element of type T to default_label, and the latter should never be used as label for an existing edge
                        AdjMatrix(int size, T default_label, vector<T>& vertex_labels):num_of_vertices(size),default_label(default_label){
                                array = new T[size*size]();
                                this->vertex_labels = std::move(vertex_labels);
                        }
                        AdjMatrix(int size, T default_label):num_of_vertices(size), default_label(default_label){
                                array = new T[size*size]();
                        }
                        ~AdjMatrix(){
                                if(num_of_vertices > -1) delete[] array;
                        }
                        AdjMatrix(){
                                num_of_vertices = -1;
								default_label = T();
                        }
                        AdjMatrix(AdjMatrix<T>&& b){
                                this->num_of_vertices = b.num_of_vertices;
                                this->default_label = b.default_label;
                                this->array = b.array;
                                b.array = nullptr;
                                this->vertex_labels = std::move(b.vertex_labels);
                        }
                        AdjMatrix<T>& operator=(AdjMatrix<T>&& b){
                                this->num_of_vertices = b.num_of_vertices;
                                this->default_label = b.default_label;
                                this->array = b.array;
                                b.array = nullptr;
                                this->vertex_labels = std::move(b.vertex_labels);
                                return *this;
                        }
                        void addEdge(int v, int u, T label, bool bothDirections=false){
                                array[v*num_of_vertices+u] = label;
                                if(bothDirections) addEdge(u,v,label);
                        }
                        void delEdge(int v, int u, bool bothDirections=false){
                                addEdge(v,u,default_label, bothDirections);
                        }
                        const T& getEdge(int v, int u) const{
                                return array[v*num_of_vertices+u];
                        }
                        const T& operator()(int v, int u) const{
                                return getEdge(v,u);
                        }
                        bool hasVertexLabels() const{
                                return vertex_labels.size()!=0;
                        }
                        const T& getVertexLabel(int i) const{
                                if(!hasVertexLabels() || i >= num_of_vertices) throw out_of_range("No vertex labels");
                                return vertex_labels[i];
                        }
                private:
                        int num_of_vertices;
                        T default_label;
                        T* array;
                        vector<T> vertex_labels;
        };
        class AtomicType{
                public:
                        class VectorView{
                                public:
                                        
                                        VectorView(const AdjMatrix<int>& adjMatrix, const vector<int>& v): vertex_subset(v), adjMatrix(adjMatrix){}
                                        int operator[](size_t i) const{
                                                auto s = vertex_subset.size();
                                                auto row = i/s;
                                                auto col = i%s;
                                                auto iV = this->vertex_subset[row];
                                                auto jV = this->vertex_subset[col];
                                                if(iV == jV){
                                                        if(adjMatrix.hasVertexLabels()) return -adjMatrix.getVertexLabel(iV);
                                                        return -1;
                                                }
                                                return adjMatrix(iV,jV);
                                        }
                                private:
                                        const vector<int>& vertex_subset;
                                        const AdjMatrix<int>& adjMatrix;
                        };
                        template<typename Container>
                        AtomicType(const Container& vertex_subset, const AdjMatrix<int>& adjMatrix) : adj_matrix(adjMatrix){
                              this->vertex_subset.reserve(vertex_subset.size());
                              this->vertex_subset = vertex_subset;
                              this->size = vertex_subset.size();
                        }
                        
                        const VectorView getSetVector() const{
                                return VectorView(adj_matrix, vertex_subset);
                        }
                private:
                        
                        int size;
                        vector<int> vertex_subset;
                        const AdjMatrix<int>& adj_matrix;
        };
        template<class T>
        class Wrapper{
                private:
                        T& c;
                public:
                        Wrapper(T& c):c(c){}
                        void clearSetVector(){
                                c.clearSetVector();
                        }
                        const vector<int>& getSetVector() const{
                                return c.getSetVector();
                        }
        };
        class Coloring{
                public:
                        
                        Coloring(){
                                previous_color = -1;
                        }
                        Coloring(int pc): previous_color(pc){};
                        void update_previous_color(int newColor){
                                clearSetVector();
                                previous_color = newColor;
                        }
                        void add_multiset_color(int color, int size = -1){
                                if(size > 0) set_coloring.reserve(size);
                                set_coloring.push_back(color);
                        }
                        const vector<int>& getSetVector() const{
                                return set_coloring;
                        }
                        void sortSetVector(){
                                sort(set_coloring.begin(), set_coloring.end());
                        }
                        void clearSetVector(){
                                set_coloring = vector<int>();
                        }
                        int getColor() const{
                                return previous_color;
                        }
                private:
                        int previous_color;
                        vector<int> set_coloring;
        };
        
        void innerGenerateTupleMap(int i, int k, int n, vector<int>& currentTuple, TupleMap& tupleMap, InverseTupleMap& inverseTupleMap){
                if(i == k){
                        auto t = Tuple<int>(currentTuple);
                        size_t idx = tupleMap.size();   
                        inverseTupleMap[idx] = tupleMap.emplace(std::move(t), idx).first;
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
                                res.push_back(std::vector<int>());
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
        pair<vector<int>, vector<bool>> orderSortedSets(const vector<T>& sortedSets, int setSize){ //First returned value is the remap vector (that is v[i]= j means the j-th element should become the i-th) while the second is a vector containing numbers > 0 at each index where a new bucket begins
                using BucketTuple = std::tuple<int,int,int>;
                vector<int> orderMap(sortedSets.size());
                vector<bool> buckets(sortedSets.size());
                iota<int>(orderMap.begin(), orderMap.end(), 0);
                buckets[0] = 1;
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
                                buckets[tupleB+b] = true;
                                //std::cout << "------- SottoBucket: (" << tupleB+b << ", " << tupleB+e << ", " << sortingIndex+1 << ")" << std::endl;
                                if(e - b > 1 && sortingIndex + 1 != setSize) bucketsToSort.push({tupleB+b, tupleB+e, sortingIndex+1});
                                b = e;
                        }
                }
                return {orderMap, buckets};
        }
        int updateColoring(const vector<int>& remap, const vector<bool>& buckets, vector<Coloring>& colorings){
                int k = -1, s = remap.size();
                int cnt = 0;
                for(int i = 0; i < s; i++){
                        if(buckets[i]){
                                ++cnt; 
                                ++k;
                        }
                        colorings[remap[i]].update_previous_color(k);
                }
                return cnt;
        }
        //Returns the number of new colors
        int updateColoring(const vector<int>& remap, const vector<int>& buckets, vector<Coloring>& colorings, queue<int>& colorQueue, int lastColor = -1){
                int k = lastColor, s = remap.size();
                int size = 0;
                int maxSize = -1, maxColor = -1;
                vector<int> colors;
                for(int i = 0; i < s; i++){
                        if(buckets[i] > 0){
                                if(colors.size() > 0){
                                        colors.push_back(k);
                                        if(size > maxSize){
                                                maxSize = size;
                                                maxColor = k;
                                        }
                                }
                                size = 0; 
                                ++k;
                        }
                        colorings[remap[i]].update_previous_color(k);
                        ++size;
                }
                for(const auto& color: colors){
                        if(color != maxColor){
                                colorQueue.push(color);
                        }
                }
                return colors.size();
        }
        
        
        Coloring& prepareElementColoring(const vector<Coloring>& firstColoring, TupleMap& tm, const InverseTupleMap& itm, vector<Coloring>& tuple_coloring, int i, int n, int k){
                auto& color = tuple_coloring[i];
                const auto& tuple = itm[i]->first;
                vector<vector<int>> tmpMultiset;
                tmpMultiset.reserve(n);
                for(int v = 0; v < n; v++){
                        tmpMultiset.push_back(vector<int>(k));
                        for(int t = 0; t < k; t++){
                                auto tmpTuple = tuple.modify(t, v);
                                int tmpTupleIndex = tm[tmpTuple];
                                tmpMultiset.back()[t] = tuple_coloring[tmpTupleIndex].getColor();
                        }
                }
                sort(tmpMultiset.begin(), tmpMultiset.end());
                if(k == 1) color.add_multiset_color(firstColoring[i].getColor(), n*k);
                for(const auto& innerVector:tmpMultiset){
                        for(const auto& el: innerVector){
                                color.add_multiset_color(el, n*k);
                        }
                }
                return color;
        }
        
        void disposeElementColoring(Wrapper<Coloring>& el){
                el.clearSetVector();
        }
        
        template<class T>
        bool orderSetOfSetsBuckets(vector<int>& remap, vector<bool>& buckets, vector<Coloring>& tuple_coloring, std::function<void(Wrapper<T>&)> disposeElement, std::function<T&(int)> prepareElement, int setSize){
                bool finished = true;
                int s = buckets.size();
                vector<Wrapper<T>> tmp;
                int last_i = 0;
                int i = 0;
                while(i < s){
                        do{
                                T& el = prepareElement(remap[i]);
                                tmp.push_back(Wrapper<T>(el));
                                ++i;
                        }while(i < s && !buckets[i]);
                        auto new_res = orderSortedSets(tmp, setSize);
                        auto& new_remap = new_res.first;
                        auto& new_buckets = new_res.second;
                        int bucketCounter = 0;
                        auto new_remap_size = new_remap.size();
                        for(size_t j = 0; j < new_remap_size; j++){
                                new_remap[j] = remap[last_i + new_remap[j]];
                        }
                        for(auto& el: tmp) disposeElement(el);
                        for(int j = last_i; j < i; j++){
                                remap[j] = new_remap[j-last_i];
                                if(new_buckets[j-last_i]){
                                        bucketCounter++;
                                        buckets[j] = true;
                                }
                                if(bucketCounter > 1) finished = false;
                        }
                        last_i = i;
                        tmp.clear();
                }
                updateColoring(remap, buckets, tuple_coloring);
                return finished;
        }
        //Edge labels should start from 1
        TupleMap k_WL(const std::vector<GraphNode>& v, int k, bool hasVertexLabels){
                int n = v.size();
                TupleMap tm;
                InverseTupleMap itm;
                vector<Coloring> tuple_coloring((int)pow(n,k));
                vector<Coloring> firstColoring;
                pair<vector<int>,vector<bool>> res;
                
                {       //The block is used to send adj_matrix and atps out of scope when not needed anymore, so as to free some memory
                        AdjMatrix<int> adj_matrix(n, 0);
                        if(hasVertexLabels){
                                vector<int> vertex_labels;
                                vertex_labels.reserve(n);
                                for(const auto& el: v){
                                         vertex_labels.push_back(el.color);
                                }
                                adj_matrix = std::move(AdjMatrix<int>(n,0,vertex_labels));
                        }
                        for(const auto& el: v){
                                int vIdx = el.idx;
                                for(const auto& adj: el.adj_list){
                                        adj_matrix.addEdge(vIdx, adj.first, adj.second);
                                }
                        }
                
                        generateTupleMap(n, k, tm, itm);
                
                        auto numberOfTuples = itm.size();
                        vector<AtomicType> atp;
                        atp.reserve(numberOfTuples);
                        for(const auto& el: itm){
                                atp.push_back(AtomicType(el->first, adj_matrix));
                        }
                        res = orderSortedSets(atp, k*k);
                        auto& atp_remap = res.first;
                        auto& atp_buckets = res.second;
                        
                        //Now the idea is going through the atps in remapping order,
                        //and initialize the color of each tuple based on the index of the bucket the tuple is in
                        
        
                        updateColoring(atp_remap, atp_buckets, tuple_coloring);
                
                
                        if(k == 1) firstColoring = tuple_coloring;
                }
                
                
                //After this, the main part of the algorithm, comprised of computing the coloring of adjacents, ordering the
                //tuples by old_coloring^multiset and then update their old coloring.
                //One should stop when the orbits are the same before and after a round. This last part is gonna be tricky for sure
                
                bool finished = false;
                while(!finished){
                        auto& remap = res.first;
                        auto& buckets = res.second;
                        finished = orderSetOfSetsBuckets<Coloring>(remap, buckets, tuple_coloring, std::function<void(Wrapper<Coloring>&)>(disposeElementColoring), std::function<Coloring&(int)>([&firstColoring, &tm, &itm, &tuple_coloring, n, k](int i)->Coloring&{return prepareElementColoring(firstColoring, tm, itm, tuple_coloring, i, n, k);}), n*k);
                }
                for(auto& el: tm){
                        el.second = tuple_coloring[el.second].getColor();
                }
                return tm;
        }

}
