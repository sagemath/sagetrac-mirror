/*
 The MIT License (MIT)

 Copyright (c) 2015 Yuki Kawata

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 the Software, and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __GRAPH_DIAMETER_H__
#define __GRAPH_DIAMETER_H__

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <stack>
#include <algorithm>
#include <sys/time.h>

class GraphDiameter {
    public:
    
    GraphDiameter() : V(0), diameter(0), numBFS(0), time(0) {}
    ~GraphDiameter() {}
    int GetDiameter(const std::vector <std::pair<int, int> > &edges, int num_double_sweep = numDefaultDoubleSweep);
    int GetDiameter(const char *filename, int num_double_sweep = numDefaultDoubleSweep);
    void PrintStatistics(void);
    
    private:
    
    static const int numDefaultDoubleSweep = 10;
    int V;
    int diameter;
    int numBFS;
    double time;
    
    double GetTime(void) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
    
    int GetRandom(void) {
        static unsigned x = 123456789;
        static unsigned y = 362436039;
        static unsigned z = 521288629;
        static unsigned w = 88675123;
        unsigned t;
        
        t = x ^ (x << 11);
        x = y;
        y = z;
        z = w;
        w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
        
        return w % V;
    }
};

int GraphDiameter::GetDiameter(const std::vector <std::pair<int, int> > &edges, int num_double_sweep) {
    // Prepare the adjacency list
    std::vector <std::vector <int> > graph;
    std::vector <std::vector <int> > rgraph;
    {
        for (size_t i = 0; i < edges.size(); i++) {
            int from = edges[i].first;
            int to = edges[i].second;
            
            V = std::max(V, from + 1);
            V = std::max(V, to + 1);
        }
        
        graph.resize(V);
        rgraph.resize(V);
        
        for (size_t i = 0; i < edges.size(); i++) {
            int from = edges[i].first;
            int to = edges[i].second;
            
            graph[from].push_back(to);
            rgraph[to].push_back(from);
        }
    }
    
    // Decompose the graph into strongly connected components
    time = -GetTime();
    std::vector <int> scc(V);
    {
        int num_visit = 0, num_scc = 0;
        std::vector <int> ord(V, -1);
        std::vector <int> low(V);
        std::vector <bool> in(V, false);
        std::stack <int> s;
        std::stack <std::pair<int, int> > dfs;
        
        for (int i = 0; i < V; i++) {
            if (ord[i] != -1) continue;
            
            dfs.push(std::make_pair(i, -1));
            
            while (!dfs.empty()) {
                int v = dfs.top().first;
                int index = dfs.top().second;
                
                dfs.pop();
                
                if (index == -1) {
                    ord[v] = low[v] = num_visit++;
                    s.push(v);
                    in[v] = true;
                } else {
                    low[v] = std::min(low[v], low[graph[v][index]]);
                }
                
                for (index++; index < (int)graph[v].size(); index++) {
                    int w = graph[v][index];
                    
                    if (ord[w] == -1) {
                        dfs.push(std::make_pair(v, index));
                        dfs.push(std::make_pair(w, -1));
                        break;
                    } else if (in[w] == true) {
                        low[v] = std::min(low[v], ord[w]);
                    }
                }
                
                if (index == (int)graph[v].size() && low[v] == ord[v]) {
                    while (true) {
                        int w = s.top();
                        
                        s.pop();
                        in[w] = false;
                        scc[w] = num_scc;
                        
                        if (v == w) break;
                    }
                    
                    num_scc++;
                }
            }
        }
    }
    
    // Compute the diameter lower bound by the double sweep algorithm
    int qs, qt;
    std::vector <int> dist(V, -1);
    std::vector <int> queue(V);
    {
        for (int i = 0; i < num_double_sweep; i++) {
            int start = GetRandom();
            
            // forward BFS
            qs = qt = 0;
            dist[start] = 0;
            queue[qt++] = start;
            
            while (qs < qt) {
                int v = queue[qs++];
                
                for (size_t j = 0; j < graph[v].size(); j++) {
                    if (dist[graph[v][j]] < 0) {
                        dist[graph[v][j]] = dist[v] + 1;
                        queue[qt++] = graph[v][j];
                    }
                }
            }
            
            for (int j = 0; j < qt; j++) dist[queue[j]] = -1;
            
            // barkward BFS
            start = queue[qt - 1];
            qs = qt = 0;
            dist[start] = 0;
            queue[qt++] = start;
            
            while (qs < qt) {
                int v = queue[qs++];
                
                for (size_t j = 0; j < rgraph[v].size(); j++) {
                    if (dist[rgraph[v][j]] < 0) {
                        dist[rgraph[v][j]] = dist[v] + 1;
                        queue[qt++] = rgraph[v][j];
                    }
                }
            }
            
            diameter = std::max(diameter, dist[queue[qt - 1]]);
            
            for (int j = 0; j < qt; j++) dist[queue[j]] = -1;
        }
    }
    
    // Order vertices
    std::vector <std::pair<long long, int> > order(V);
    {
        for (int v = 0; v < V; v++) {
            int in = 0, out = 0;
            
            for (size_t i = 0; i < rgraph[v].size(); i++) {
                if (scc[rgraph[v][i]] == scc[v]) in++;
            }
            
            for (size_t i = 0; i < graph[v].size(); i++) {
                if (scc[graph[v][i]] == scc[v]) out++;
            }
            
            // SCC : reverse topological order
            // inside an SCC : decreasing order of the product of the indegree and outdegree for vertices in the same SCC
            order[v] = std::make_pair(((long long)scc[v] << 32) - in * out, v);
        }
        
        std::sort(order.begin(), order.end());
    }
    
    // Examine every vertex
    std::vector <int> ecc(V, V);
    {
        for (int i = 0; i < V; i++) {
            int u = order[i].second;
            
            if (ecc[u] <= diameter) continue;
            
            // Refine the eccentricity upper bound
            int ub = 0;
            std::vector <std::pair<int, int> > neighbors;
            
            for (size_t j = 0; j < graph[u].size(); j++) neighbors.push_back(std::make_pair(scc[graph[u][j]], ecc[graph[u][j]] + 1));
            
            sort(neighbors.begin(), neighbors.end());
            
            for (size_t j = 0; j < neighbors.size(); ) {
                int component = neighbors[j].first;
                int lb = V;
                
                for (; j < neighbors.size(); j++) {
                    if (neighbors[j].first != component) break;
                    lb = std::min(lb, neighbors[j].second);
                }
                
                ub = std::max(ub, lb);
                
                if (ub > diameter) break;
            }
            
            if (ub <= diameter) {
                ecc[u] = ub;
                continue;
            }
            
            // Conduct a BFS and update bounds
            numBFS++;
            qs = qt = 0;
            dist[u] = 0;
            queue[qt++] = u;
            
            while (qs < qt) {
                int v = queue[qs++];
                
                for (size_t j = 0; j < graph[v].size(); j++) {
                    if (dist[graph[v][j]] < 0) {
                        dist[graph[v][j]] = dist[v] + 1;
                        queue[qt++] = graph[v][j];
                    }
                }
            }
            
            ecc[u] = dist[queue[qt - 1]];
            diameter = std::max(diameter, ecc[u]);
            
            for (int j = 0; j < qt; j++) dist[queue[j]] = -1;
            
            qs = qt = 0;
            dist[u] = 0;
            queue[qt++] = u;
            
            while (qs < qt) {
                int v = queue[qs++];
                
                ecc[v] = std::min(ecc[v], dist[v] + ecc[u]);
                
                for (size_t j = 0; j < rgraph[v].size(); j++) {
                    // only inside an SCC
                    if (dist[rgraph[v][j]] < 0 && scc[rgraph[v][j]] == scc[u]) {
                        dist[rgraph[v][j]] = dist[v] + 1;
                        queue[qt++] = rgraph[v][j];
                    }
                }
            }
            
            for (int j = 0; j < qt; j++) dist[queue[j]] = -1;
        }
    }
    
    time += GetTime();
    
    return diameter;
}

int GraphDiameter::GetDiameter(const char *filename, int num_double_sweep) {
    FILE *in = fopen(filename, "r");
    
    if (in == NULL) {
        fprintf(stderr, "Can't open the graph file\n");
        return -1;
    }
    
    std::vector <std::pair<int, int> > edges;
    
    for (int from, to; fscanf(in, "%d %d", &from, &to) != EOF; ) edges.push_back(std::make_pair(from, to));
    
    fclose(in);
    
    return GetDiameter(edges, num_double_sweep);
}

void GraphDiameter::PrintStatistics(void) {
    printf("Diameter : %d\n", diameter);
    printf("#BFS : %d -> %d\n", V, numBFS);
    printf("Time : %lf sec\n", time);
}
#endif
