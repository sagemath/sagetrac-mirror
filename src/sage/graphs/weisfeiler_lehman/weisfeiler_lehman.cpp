#include "weisfeiler_lehman.h"
#include "Tuple.h"
#include <vector>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <cstring>
#include <exception>
#include <cassert>
namespace wl{
        template<typename T>
        class AdjMatrix{
                public:
                        int* getCanonicalOrdering(int* vertices, int s, int n ,int iS, int jS) const{
                            int* permutation = new int[s];
                            memcpy(permutation, vertices, s*sizeof(int));
                            int* current_best = new int[s];
                            memcpy(current_best, vertices, s*sizeof(int));
                            char* m = new char[n*n]();
                            getCanonicalOrdering(vertices, current_best, 0, s, n, m, iS, jS);
                            delete[] m;
                            delete[] permutation;
                            return current_best;
                        }
                        class VectorView{
                                public:
                                        VectorView(const AdjMatrix<int>& adjMatrix, int* v, int s, pair<int,int> special_edge): adjMatrix(adjMatrix){
                                            vertex_subset = v;
                                            this->s = s;
                                            this->special_edge = special_edge;
                                            rehash();
                                        }
                                        ~VectorView(){
                                            delete[] vertex_subset;
                                        }
                                        VectorView(const VectorView& b): adjMatrix(b.adjMatrix){
                                            vertex_subset = new int[b.s];
                                            memcpy(vertex_subset, b.vertex_subset, sizeof(int)*b.s);
                                            s = b.s;
                                            special_edge = b.special_edge;
                                            rehash();
                                        }
                                        int operator[](size_t i) const{
                                                auto row = i/s;
												auto col = i%s;
												return (*this)(row, col);
                                        }
                                        int operator()(int row, int col, bool trueIndex=false) const{
                                            return adjMatrix.getValue(row, col, this->vertex_subset, special_edge.first, special_edge.second, trueIndex);
                                        }
                                        size_t getHash() const{
                                            return hash_value;
                                        }
                                        void rehash(){
                                            hash_value = 0;
                                            int n = s*s;
                                            for(int i = 0; i < n; i++){
                                                hash_combine<int>(hash_value, (*this)[i]);
                                            }
                                        }
                                        bool operator==(const VectorView& b) const{
                                            if(s != b.s) return false;
                                            if(b.hash_value != hash_value) return false;
                                            int n = s*s;
                                            for(int i = 0; i < n;i++){
                                                if((*this)[i] != b[i]) return false;
                                            }
                                            return true;
                                        }
                                        bool operator<(const VectorView& b) const{
                                            if(s != b.s) throw std::invalid_argument("Different vertex subset size");
                                            auto n = s*s;
                                            for(int i = 0; i < n;i++){
                                                int aV = (*this)[i];
                                                int bV = b[i];
                                                if(aV < bV) return true;
                                                if(bV < aV) return false;
                                            }
                                            return false;
                                        }
                                        void printVectorView() const{
                                            for(int i = 0; i < s; i++){
                                                for(int j = 0; j < s; j++){
                                                    cout << (*this)(i,j) << " ";
                                                }
                                                cout << endl;
                                            }
                                        }
                                        void newSpecialEdge(int i, int j){
                                            special_edge.first = i;
                                            special_edge.second = j;
                                            rehash();
                                        }
                                        int* vertex_subset;
                                        int s;
                                private:
                                        size_t hash_value = 0;
                                        
                                        
                                        const AdjMatrix<int>& adjMatrix;
                                        pair<int,int> special_edge;
                        };
                        //The type T default constructor should set an element of type T to default_label, and the latter should never be used as label for an existing edge
                        AdjMatrix(int size):num_of_vertices(size){
                                array = new T[size*size]();
                        }
                        ~AdjMatrix(){
                                if(num_of_vertices > -1) delete[] array;
                        }
                        AdjMatrix(){
                                num_of_vertices = -1;
                        }
                        AdjMatrix(AdjMatrix<T>&& b){
                                this->num_of_vertices = b.num_of_vertices;
                                this->array = b.array;
                                b.array = nullptr;
                                this->minValue = b.minValue;
                        }
                        AdjMatrix<T>& operator=(AdjMatrix<T>&& b){
                                this->num_of_vertices = b.num_of_vertices;
                                this->array = b.array;
                                b.array = nullptr;
                                this->minValue = b.minValue;
                                return *this;
                        }
                        void addEdge(int v, int u, T label, bool bothDirections=false){
                                array[v*num_of_vertices+u] = label;
                                if(label < minValue) minValue = label;
                                if(bothDirections) addEdge(u,v,label);
                        }
                        const T& getEdge(int v, int u) const{
                            return array[v*num_of_vertices+u];
                        }
                        const T& operator()(int v, int u) const{
                                return getEdge(v,u);
                        }
                        int size() const{
                            return num_of_vertices;
                        }
                        template<typename U>
                        static void swap(AdjMatrix<U>& a, AdjMatrix<U>& b){
                            int t = a.num_of_vertices;
                            a.num_of_vertices = b.num_of_vertices;
                            b.num_of_vertices = t;
                            U* tempArr = a.array;
                            a.array = b.array;
                            b.array = tempArr;
                            U tempMin = a.minValue;
                            a.minValue = b.minValue;
                            b.minValue = tempMin;
                        }
                private:
                        int num_of_vertices;
                        T* array;
                        T minValue = T();
                        int compareSubgraphs(int* a, int* b, size_t s, size_t determinationPoint, int iS, int jS) const{
                            if(memcmp(a, b, s*sizeof(int)) == 0) return 0;
                            for(size_t i = 0; i < s; i++){
                                for(size_t j = 0; j < s; j++){
                                    int aV = getValue(i, j, a, iS, jS);
                                    int bV = getValue(i, j, b, iS, jS);
                                    if(aV < bV) return -1;
                                    if(bV < aV){
                                        if(i <= determinationPoint && j <= determinationPoint){
                                            if(i == determinationPoint || j == determinationPoint) return 3;
                                            else return 2;
                                        } else return 1;
                                    }
                                }
                            }
                            return 0;
                        }
                        bool areVerticesIsomorphic(int a, int b, int* v, int s, int iS, int jS) const{
                            if(v[a] == v[b]) return true;
                            if(getValue(a,a,v,iS,jS) != getValue(b,b,v,iS,jS) || getValue(a,b,v,iS,jS) != getValue(b,a,v,iS,jS)) return false;
                            for(int i = 0; i < s; i++){
                                if(i == b || i == a) continue;
                                int aV = getValue(a, i, v,iS,jS);
                                int bV = getValue(b, i, v, iS, jS);
                                if(aV != bV) return false;
                                aV = getValue(i, a, v, iS, jS);
                                bV = getValue(i, b, v, iS, jS);
                                if(aV != bV) return false;
                            }
                            return true;
                        }
                        /* OLD IMPLEMENTATION, SHOULD BE A LOT WORSE BUT IT'S NOT, NEED TO CHECK WHY
                           void getCanonicalOrdering(int* permutation, int* current_best, size_t offset, size_t s, size_t n, char* isomorphicVertices) const{
                            int c = compareSubgraphs(permutation, current_best, s, offset);
                            if(c == 2) return;
                            if(c == -1) memcpy(current_best, permutation, s*sizeof(int));
                            for(size_t i = (offset == s-2)?offset+1:offset; i < s; i++){
                                if(i != offset){
                                    int idx = permutation[i]*n+permutation[offset];
                                    char v = isomorphicVertices[idx];
                                    if(v == 0){
                                        bool r = areVerticesIsomorphic(i, offset, permutation, s);
                                        isomorphicVertices[idx] = r?1:2;
                                        if (r) continue;
                                    }else{
                                        if (v == 1) continue;
                                    }
                                }
                                std::iter_swap(permutation+i, permutation + offset);
                                getCanonicalOrdering(permutation, current_best, offset+1, s, n, isomorphicVertices);
                                std::iter_swap(permutation + i, permutation + offset);
                            }
                        }*/
                        
                        void getCanonicalOrdering(int* permutation, int* current_best, size_t offset, size_t s, size_t n, char* isomorphicVertices, int iS, int jS) const{
                            for(size_t i = (offset == s-2)?offset+1:offset; i < s; i++){
                                if(i != offset){
                                    int idx = permutation[i]*n+permutation[offset];
                                    char v = isomorphicVertices[idx];
                                    if(v == 0){
                                        bool r = areVerticesIsomorphic(i, offset, permutation, s, iS, jS);
                                        isomorphicVertices[idx] = r?1:2;
                                        if (r) continue;
                                    }else{
                                        if (v == 1) continue;
                                    }
                                }
                                std::iter_swap(permutation+i, permutation + offset);
                                int c = 0;
                                if(i != offset) c = compareSubgraphs(permutation, current_best, s, offset, iS, jS);
                                if(c == -1) memcpy(current_best, permutation, s*sizeof(int));
                                if(c != 2 && c != 3) getCanonicalOrdering(permutation, current_best, offset+1, s, n, isomorphicVertices, iS, jS);
                                std::iter_swap(permutation + i, permutation + offset);
                                if(c == 2) return;
                                if(c == 3) continue;                                
                            }
                        }
                        int getValue(int row, int col, int* vertices, int special_row, int special_col, bool trueIndex=false) const{
                            int iV = row; 
                            int jV = col;
                            if(!trueIndex){
                                iV = vertices[row];
                                jV = vertices[col];
                            }
                            if(iV == special_row && jV == special_col) return minValue-1;
                            return (*this)(iV,jV);
                        }
							
        };
        

        struct VectorView_Hash{
            size_t operator() (const AdjMatrix<int>::VectorView& vv) const {
                return vv.getHash();
            }
        };
       
        //Edge labels should start from 1
        AdjMatrix<int> populateAdjMatrix(const std::vector<GraphNode>& v, int n, bool hasVertexLabels){
            AdjMatrix<int> adj_matrix(n);
            for(const auto& el: v){
                    int vIdx = el.idx;
                    for(const auto& adj: el.adj_list){
                            adj_matrix.addEdge(vIdx, adj.first, adj.second);
                    }
            }
            for(const auto& el: v){
                     adj_matrix.addEdge(el.idx, el.idx, hasVertexLabels?-(el.color):-1);
            }
            
            return adj_matrix;
        }
        
        
        using ColorClass = std::unordered_set<std::pair<int,int>, IntPair_Hash>;
        using FingerprintMap = std::unordered_map<AdjMatrix<int>::VectorView, int, VectorView_Hash>;
        int initColorClasses(queue<ColorClass>& cc, const AdjMatrix<int>& am){
            int n = am.size();
            unordered_map<int, ColorClass> temp;
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    temp[am(i,j)].insert({i,j});
                }
            }
            int c = temp.size();
            for(auto& color_class: temp){
                cc.push(std::move(color_class.second));
            }
            return c;
        }
        void initFingerprint(FingerprintMap& fingerprint, const AdjMatrix<int>& am, int maxVertex, int offset, int limit, int* currentTuple, unordered_set<int>& used_vertices){
            int n = am.size();
            if(offset == limit){ //If a subgraph is completed
                for(const auto& i: used_vertices){
                    for(const auto& j: used_vertices){
                        AdjMatrix<int>::VectorView subgraphView(am, am.getCanonicalOrdering(currentTuple, limit, n, i,j), limit, {i,j});
                        fingerprint[subgraphView] = 0;
                    }
                }
            }else{            
                for(int k = maxVertex; k < n; k++){
                    currentTuple[offset] = k;
                    bool inserted = used_vertices.insert(k).second;
                    initFingerprint(fingerprint, am, k, offset+1, limit, currentTuple, used_vertices);
                    if(inserted) used_vertices.erase(k);
                }
            }
        }
        void createFingerprint(FingerprintMap& fingerprint, const AdjMatrix<int>& am, int i, int j, int maxVertex, int offset, int limit, int* currentTuple, bool usedI=false, bool usedJ=false){
            if(offset > limit-2 && !usedI && !usedJ) return;
            if((maxVertex > i && !usedI) || (maxVertex > j && !usedJ)) return;
            int n = am.size();
            if(offset == limit){ //If a subgraph is completed
                if(!usedI || !usedJ) return;
                AdjMatrix<int>::VectorView subgraphView(am, am.getCanonicalOrdering(currentTuple, limit, n, i,j), limit, {i,j});
                fingerprint[subgraphView]++;
                if(i > 8000) subgraphView.printVectorView();
            }else{            
                for(int k = maxVertex; k < n; k++){
                    currentTuple[offset] = k;
                    createFingerprint(fingerprint, am, i, j, k, offset+1, limit, currentTuple, k == i?true:usedI, (k == j && (i!=j || usedI))?true:usedJ);
                }
            }
            
        }
        int clearFingerprint(FingerprintMap& fingerprint){
            for(auto& v: fingerprint){
                v.second = 0;
            }
            return fingerprint.size();
        }
        void printFingerprint(FingerprintMap& fingerprint){
            for(const auto& el: fingerprint){
                cout << "Subset: ";
                for(int i = 0; i < el.first.s; i++){
                    cout << el.first.vertex_subset[i] << " ";
                }
                cout << endl;
                el.first.printVectorView();
                cout << endl << endl;
            }
        }
        
        unordered_map<int, vector<pair<int,int>>> k_WL(const std::vector<GraphNode>& v, int k, bool hasVertexLabels){
            //for(const auto& el: v){
                //cout << "IDX: " << el.idx << " || COLOR: " << el.color << endl;
                //cout << "-----[";
                //for(const auto& el2: el.adj_list){
                    //cout << "(" << el2.first << "," << el2.second << "), ";
                //}
                //cout << "]" << endl;
            //}
            int n = v.size();
            AdjMatrix<int> adjMatrix = populateAdjMatrix(v, n, hasVertexLabels);
            AdjMatrix<int> new_adjMatrix = AdjMatrix<int>(n);
            queue<ColorClass> color_classes;
            initColorClasses(color_classes, adjMatrix);
            queue<ColorClass> new_color_classes;
            unordered_set<int> used_vertices;
            bool finished = false;
            while(!finished){
                    finished = true;
                    FingerprintMap fingerprint;
                    int* tempVector = new int[k+1];
                    unordered_map<Tuple<int>, ColorClass> fingerprintsDB;
                    initFingerprint(fingerprint, adjMatrix, 0, 0, k+1, tempVector, used_vertices);
                    int c = 0;
                    while(!color_classes.empty()){
                        auto cc = color_classes.front();
                        color_classes.pop();
                        fingerprintsDB.clear();
                        //printFingerprint(fingerprint);
                        
                        for(const auto& edge: cc){
                            int s = clearFingerprint(fingerprint);
                            createFingerprint(fingerprint, adjMatrix, edge.first, edge.second, 0, 0, k+1, tempVector);
                            assert(s == (int)fingerprint.size());
                            Tuple<int> t(s);
                            int c = 0;
                            for(const auto& el: fingerprint){
                                t[c++] = el.second;
                            }
                            t.rehash();
                            fingerprintsDB[t].insert(edge);
                        }
                        int oldc = c;
                        for(auto& color_class: fingerprintsDB){
                            for(const auto& edge: color_class.second){
                                new_adjMatrix.addEdge(edge.first, edge.second, c);
                            }
                            c++;
                            new_color_classes.push(std::move(color_class.second));
                        }
                        if(c-oldc > 1) finished = false;
                        
                    }
                    
                    delete[] tempVector;    
                    AdjMatrix<int>::swap(adjMatrix, new_adjMatrix);
                    std::swap(color_classes, new_color_classes);                        
            }
            unordered_map<int, vector<pair<int,int>>> result;
            int c = 0;
            while(!color_classes.empty()){
                auto v = std::move(color_classes.front());
                color_classes.pop();
                result[c++] = vector<pair<int,int>>(v.begin(), v.end());
            }
            //for(const auto& el: result){
                //cout << el.first << ":" << endl;
                //for(const auto& el2: el.second){
                    //cout << "    (" << el2.first << ", " << el2.second << ")" << endl;
                //}
            //}
            return result;
        }

}
