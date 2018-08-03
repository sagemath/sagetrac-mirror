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
#include <chrono>
#include <exception>
#include <cassert>
namespace wl{
        template<typename T>
        class AdjMatrix{
                public:
                        class VectorView{
                                /*
                                 * This class allows to define a subset of vertices of a graph represented by an AdjMatrix, and create a view
                                 * of the AdjMatrix representing the induced subgraph's adjacency matrix.
                                 * Notice that no new adjacency matrix is actually created, this class is just a view that, when iterated, gives
                                 * the impression of having an adjacency matrix for the subgraph.
                                 * It is also possible to define a "special edge", that is an ordered pair of vertices, an entry in the subgraph's 
                                 * adjacency matrix, that will always be individualised, specifically will always have a unique value, whatever its
                                 * actual value in the main AdjMatrix is.
                                 * An hash is precalculated for the view at the moment of creation/copy/move, this means that the rehash() method
                                 * needs to be called if the underlying adjacency matrix changes and the view needs to remain valid.
                                 */
                                public:
                                        VectorView(const AdjMatrix<int>& adjMatrix, int* v, int s, pair<int,int> special_edge): adjMatrix(adjMatrix){
                                            vertex_subset = v;
                                            this->s = s;
                                            this->special_edge = std::move(special_edge);
                                            rehash();
                                        }
                                        VectorView(const AdjMatrix<int>& adjMatrix, int* v, int s): adjMatrix(adjMatrix){
                                            vertex_subset = v;
                                            this->s = s;
                                            this->special_edge = {-1,-1};
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
                                            hash_value = b.hash_value;
                                        }
                                        VectorView(VectorView&& b): adjMatrix(b.adjMatrix){
                                            vertex_subset = b.vertex_subset;
                                            b.vertex_subset = nullptr;
                                            s = b.s;
                                            special_edge = std::move(b.special_edge);
                                            hash_value = b.hash_value;
                                        }
                                        //Sequential access, as if the adjacency matrix was stored in a one-dimensional array
                                        int operator[](size_t i) const{
                                                auto row = i/s;
												auto col = i%s;
												return (*this)(row, col);
                                        }
                                        /*
                                         * Matrix access, allows to specify the row and column of the entry one wants to retrieve.
                                         * Notice that the indices provided are the indices of the entry in the subgraph's adjacency matrix, not in
                                         * the main matrix.
                                         * E.g: if the subgraph has the following vertex set (0,5,2,3,6) retrieving the entry with row = 1 and col=3
                                         * will return the entry for (5,3).
                                         * If one wants to actually retrieve the entry for (row, col), set the parameter trueIndex to true
                                         */
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
                                        
                                private:
                                        size_t hash_value = 0;
                                        int* vertex_subset;
                                        int s;
                                        const AdjMatrix<int>& adjMatrix;
                                        pair<int,int> special_edge;
                        };
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
                            AdjMatrix<U> t(std::move(a));
                            a = std::move(b);
                            b = std::move(t);
                        }
                        
                        /*
                         * Given a subset of vertices of "this", possibly with repetitions, the method returns a vertex ordering that
                         * lexicographically minimizes the induced subgraph's adjacency matrix, taking into account the eventual individualisation
                         * of a pair (iS, jS) as described in the comments to the VectorView class
                         */
                        int* getCanonicalOrdering(int* vertices, int s, int n ,int iS, int jS, char* m) const{
                            int* current_best = new int[s];
                            memcpy(current_best, vertices, s*sizeof(int));
                            getCanonicalOrdering(vertices, current_best, 0, s, n, m, iS, jS);
                            return current_best;
                        }
                private:
                        int num_of_vertices;
                        T* array;
                        T minValue = T(); //Value with which to individualise entries, for VectorView's usage
                        
                        /*
                         * Private method used to lexicographically compare the adjacency matrices of the subgraphs induced by
                         * two ordered (multi)sets of vertices of "this".
                         * It takes into account the individualisation of an entry as described in the comments to the VectorView class.
                         * It also accepts a determination point, that is an index smaller than "s" that tells the algorithm the index
                         * up to which the ordering is "fixed", and after which it can still change without requiring the calling method to backtrack.
                         * 
                         * It returns 0 if the two adjacency matrices are equal.
                         * If the adjacency matrix induced by "a" is lexicographically smaller than the one induced by "b" the method will return a negative number:
                         *      > "-1" if the first value by which the two matrices differ has both the row and column indices larger than the determination point
                         *      > "-2" if the first value by which the two matrices differ has both the row and column indices smaller than the determination point
                         *      > "-3" otherwise
                         * The same rules about the modulo returned apply if the adjacency matrix of "b" is the lexicographically larger, but the value returned will be positive
                         */
                        int compareSubgraphs(int* a, int* b, size_t s, size_t determinationPoint, int iS, int jS) const{
                            if(memcmp(a, b, s*sizeof(int)) == 0) return 0; //Check if the sets of vertices are equal, in that case the adjacency matrices will be too, so we can return 0
                            for(size_t i = 0; i < s; i++){
                                for(size_t j = 0; j < s; j++){
                                    int aV = getValue(i, j, a, iS, jS);
                                    int bV = getValue(i, j, b, iS, jS);
                                    if(aV < bV){
                                        if(i <= determinationPoint && j <= determinationPoint){
                                            if(i == determinationPoint || j == determinationPoint) return -3;
                                            else return -2;
                                        } else return -1;
                                    }
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
                        
                        /*
                         * This method checks if two vertices of the subgraph induced by "v" (with indices "a" and "b" in the ordered set "v")
                         * are in the same orbit of its automorphism group. It takes into account the individualisation of an entry as described in the comments
                         * to the VectorView class.
                         */
                        bool areVerticesAutomorphic(int a, int b, int* v, int s, int iS, int jS) const{
                            if(v[a] == v[b]) return true; //If they are the same vertex, they are in the same orbit
                            //Otherwise, check that they have equal entries in their columns and in their rows (the induced matrix's ones) but at the column and row's intersection
                            //And check that they have swappable values at the intersections. If all of this holds, they are in the same orbit.
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
                        
                        /*
                         * Inner recursive method using to compute the canonical ordering as described in the comments to the public method "getCanonicalOrdering"
                         */
                        void getCanonicalOrdering(int* permutation, int* current_best, size_t offset, size_t s, size_t n, char* isomorphicVertices, int iS, int jS) const{
                            //We act on the element in position "offset", which we'll call the current element
                            //For each element in a position between "offset" and the end of the vertex set, try to swap it with the current element and recurse, then swap it back and try the successive element in the set
                            if(offset == s-1) return;
                            for(size_t i = (offset == s-2)?offset+1:offset; i < s; i++){ //Usually we "swap" the current elements with itself to allow trying different permutations of the rest of the set. If the current element is the second to last one, though, there's no other permutation doable on the rest of the set, so it's useless to swap it with itself, and that's why i starts from offset+1
                                if(i != offset){ //Again, the fake "swap" of the current element with itself is always allowed to allow permutations of the rest of the set, even if "permutation[offset]" is clearly in the same orbit as "permutation[offset]"
                                    int idx = permutation[i]*n+permutation[offset];
                                    char v = isomorphicVertices[idx];
                                    if(v == 0){ //If we have no records on the orbits of these two vertices, check if they are the same and record it
                                        bool r = areVerticesAutomorphic(i, offset, permutation, s, iS, jS);
                                        isomorphicVertices[idx] = r?1:2;
                                        if (r) continue;
                                    }else{
                                        if (v == 1) continue; //Otherwise if they are in the same orbit, avoid swapping them, it would be fruitless
                                    }
                                }
                                std::iter_swap(permutation+i, permutation + offset);
                                int c = 0;
                                if(i != offset) c = compareSubgraphs(permutation, current_best, s, offset, iS, jS);
                                if(c < 0) memcpy(current_best, permutation, s*sizeof(int)); //If we found a better ordering, save it
                                if(c != 2 && c != 3) getCanonicalOrdering(permutation, current_best, offset+1, s, n, isomorphicVertices, iS, jS); //If we found that the best ordering found until now is smaller than the current ordering, due to a choice we have already made in previous recursions and won't be able to undo in this branch, don't continue recursing
                                std::iter_swap(permutation + i, permutation + offset);
                                if(c == 2) return; //If this ordering is larger due to a previous choice, backtrack
                                if(c == 3) continue; //If it's due to the current choice, the last swap we made, than simply continue so we can try another one.
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
       
        //Support method used to create the main adjacency matrix from the structure received as a parameter.
        //Vertex labels will be saved as negative
        //N.B: Edge and vertex labels should start from 1 and always be positive
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
        using Fingerprint = std::map<int, int>;
        
        /*
         * Iterate over the matrix to divide the pairs of vertices into equivalency classes based 
         * on their entry's value. These will be called color classes from now on
         */
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
        
        void innerCreateFingerprint(FingerprintMap& fingerprints, Fingerprint& fingerprint, const AdjMatrix<int>& am, int i, int j, int maxVertex, int offset, int limit, int* currentTuple, char* m){
            int n = am.size();
            if(offset == limit){ //If we have completely determined a vertex set
                //Create the corresponding view, based on the order given by "getCanonicalOrdering"
                AdjMatrix<int>::VectorView subgraphView(am, am.getCanonicalOrdering(currentTuple, limit, n, i,j, m), limit, {i,j});
                //Clear the support matrix used by "getCanonicalOrdering" to keep track of automorphisms. There's no need to zero out the entire matrix, since only the entries relative to the vertices in the currently computed vertex set can be different than 0
                for(int k1 = 0; k1 < limit; k1++){
                    for(int k2 = 0; k2 < limit; k2++){
                        m[currentTuple[k1]*am.size()+currentTuple[k2]] = 0;
                    }
                }
                auto& f = fingerprints[std::move(subgraphView)];
                if(f == 0){ //If this isomorphism class was never encountered before, the value of f will be 0, so we can give it a unique identifier (the number of current classes found will be enough)
                    f = fingerprints.size();
                    fingerprint[f] = 1; //And put the counter of that identifier in the fingerprint 1
                }else{
                    fingerprint[f]++; //Otherwise just increase the counter
                }
            }else{
                //With this we try all possible vertices. Since we are only interested in seeing any vertex set only once, without repetitions, and since the final ordering will be decided by getCanonicalOrdering, we can choose an arbitrary initial ordering
                //and stick with it to avoid repeating vertex sets. In our case, each vertex in the vertex set cannot be smaller than the previous ones.
                for(int k = maxVertex; k < n; k++){
                    currentTuple[offset] = k;
                    innerCreateFingerprint(fingerprints, fingerprint, am, i, j, k, offset+1, limit, currentTuple, m);
                }
            }
        }
        
        /*
         * Method to create the fingerprint for a pair or vertices (or edge), from now on called "e".
         * To do so every possible subgraph with a vertex set that includes the members of "e" and has cardinality "limit" is generated,
         * the canonical ordering how its vertex set computed, and then a counter for its isomorphism class incremented. All of this is done while individualising "e" as described in the comments to the VectorView class.
         * This will tell us how many times each possible isomorphism class (whose member have the above mentioned properties) "contains" our edge, and will allow us to distinguish them.
         */
        void createFingerprint(FingerprintMap& fingerprints, Fingerprint& fingerprint, const AdjMatrix<int>& am, int i, int j, int* currentTuple, int limit, char* m){
            //Since the vertex set will need to include "i" and "j", we can just put them as the first two members and then fill the rest of the spots, since the best ordering will be then found by "getCanonicalOrdering"
            currentTuple[0] = i;
            currentTuple[1] = j;
            innerCreateFingerprint(fingerprints, fingerprint, am, i, j, 0, 2, limit, currentTuple, m);
        }
        
        unordered_map<int, vector<pair<int,int>>> k_WL(const std::vector<GraphNode>& v, int k, bool hasVertexLabels){
            int n = v.size();
            //Create two copies of the adjacency matrix and of the color classes' queue, since during
            //the refining rounds one will be used as a reference, one will be written on, and then they will be swapped
            //before proceding with the next round
            AdjMatrix<int> adjMatrix = populateAdjMatrix(v, n, hasVertexLabels);
            AdjMatrix<int> new_adjMatrix = AdjMatrix<int>(n);
            queue<ColorClass> color_classes;
            queue<ColorClass> new_color_classes;
            
            unordered_set<int> used_vertices;
            map<int,int> fingerprint;
            char* automorphism_map = new char[n*n]();
            int* tempVector = new int[k+1];
            bool finished = false;
            
            initColorClasses(color_classes, adjMatrix);
            while(!finished){ //While at least one color class has been refined during the last round
                    finished = true;
                    FingerprintMap fingerprint_map; //Structure holding the subgraphs' adjacency matrices and the counters
                    unordered_map<vector<int>, ColorClass, IntVector_Hash> fingerprintsDB; //Structure associating each fingerprint to the corresponding edges
                    int c = 0;
                    while(!color_classes.empty()){ //For each color class
                        auto cc = color_classes.front();
                        color_classes.pop();
                        fingerprintsDB.clear(); //Clear the DB, since there can't be any bleeding between color classes
                        for(const auto& edge: cc){//For each pair of vertices (or edge) in the color class
                            createFingerprint(fingerprint_map, fingerprint, adjMatrix, edge.first, edge.second, tempVector, k+1, automorphism_map); //Create a fingerprint as described in the comments to the method
                            
                            //Notice how we don't need to generate all of the possible isomorphism classes beforehand, since if two vertices will be in the same color class, the first one analyzed will add all the needed
                            //isomorphism classes to "fingerprint_map", and the second one will use them. If the second one uses classes that the first one did not add, then they clearly won't be able to be in the same color class.
                            
                            //Save the fingerprint into an appropriately sized sorted vector
                            vector<int> t(2*fingerprint.size());
                            int idx = 0;
                            for(const auto& el: fingerprint){
                                    t[idx++] = el.first;
                                    t[idx++] = el.second;
                            }
                            //Add the edge to the correct equivalency class (which will later indicate its color class too)
                            fingerprintsDB[std::move(t)].insert(edge);
                            fingerprint.clear();
                        }
                        int oldc = c;
                        for(auto& color_class: fingerprintsDB){
                            //Create the new color classes based on the equivalency classes defined by the fingerprints, add them to the queue for the next round and write the index of the color class in the second main adjacency matrix
                            for(const auto& edge: color_class.second){
                                new_adjMatrix.addEdge(edge.first, edge.second, c);
                            }
                            c++;
                            new_color_classes.push(std::move(color_class.second));
                        }
                        if(c-oldc > 1) finished = false; //If more than one new class was generated by a single color class, it means the edges were split between them and the color class was refined, thus we need another refinement round
                        
                    }
                    AdjMatrix<int>::swap(adjMatrix, new_adjMatrix);
                    std::swap(color_classes, new_color_classes);                        
            }
            delete[] tempVector;  
            delete[] automorphism_map;
            unordered_map<int, vector<pair<int,int>>> result;
            int c = 0;
            //Put the result in a format apt for returning
            while(!color_classes.empty()){
                auto v = std::move(color_classes.front());
                color_classes.pop();
                result[c++] = vector<pair<int,int>>(v.begin(), v.end());
            }
            return result;
        }

}
