#include <boost/tuple/tuple.hpp>
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include "TD_combinations.hpp"
#include "TD_lower_bounds.hpp"
#include "TD_seperator_algorithm.hpp"
#include "TD_elimination_orderings.hpp"
#include "TD_misc.hpp"


#ifndef TD_STRUCT_VERTEX
#define TD_STRUCT_VERTEX

struct Vertex{
    unsigned int id;
};

#endif

typedef boost::adjacency_list<boost::setS, boost::listS, boost::undirectedS, Vertex> TD_graph_t;

struct bag{
    std::set<unsigned int> bag;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, bag> TD_tree_dec_t;

#include "sage_tdlib.hpp"


void make_tdlib_graph(TD_graph_t &G, std::vector<unsigned int> &V, std::vector<unsigned int> &E){
    unsigned int max = 0;
    for(unsigned int i = 0; i < V.size(); i++)
        max = (V[i]>max)? V[i] : max;

    std::vector<TD_graph_t::vertex_descriptor> idxMap(max+1);

    for(unsigned int i = 0; i < V.size(); i++){
        idxMap[V[i]] = boost::add_vertex(G);
        G[idxMap[V[i]]].id = V[i];
    }

    if(E.size() != 0){
        for(unsigned int j = 0; j < E.size()-1; j++){
            boost::add_edge(idxMap[E[j]], idxMap[E[j+1]], G);
            j++;
        }
    }
}

void make_tdlib_decomp(TD_tree_dec_t &T, std::vector<std::vector<int> > &V, std::vector<unsigned int> &E){
    std::vector<TD_tree_dec_t::vertex_descriptor> idxMap(V.size()+1);

    for(unsigned int i = 0; i < V.size(); i++){
        idxMap[i] = boost::add_vertex(T);
        std::set<unsigned int> bag;
        for(unsigned int j = 0; j < V[i].size(); j++)
            bag.insert((unsigned int) V[i][j]);
        T[idxMap[i]].bag = bag;
    }

    if(E.size() != 0){
        for(unsigned int j = 0; j < E.size()-1; j++){
            boost::add_edge(idxMap[E[j]], idxMap[E[j+1]], T);
            j++;
        }
    }

}

void make_sage_graph(TD_graph_t &G, std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    boost::graph_traits<TD_graph_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        V_G.push_back(G[*vIt].id);

    boost::graph_traits<TD_graph_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        E_G.push_back(G[boost::source(*eIt, G)].id);
        E_G.push_back(G[boost::target(*eIt, G)].id);
    }
}

void make_sage_decomp(TD_tree_dec_t &T, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    std::map<boost::graph_traits<TD_tree_dec_t>::vertex_descriptor, unsigned int> vertex_map;
    boost::graph_traits<TD_tree_dec_t>::vertex_iterator tIt, tEnd;
    unsigned int id = 0;
    
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        vertex_map.insert(std::pair<boost::graph_traits<TD_tree_dec_t>::vertex_descriptor, unsigned int>(*tIt, id++));
        std::vector<int> bag;
        for(std::set<unsigned int>::iterator sIt = T[*tIt].bag.begin(); sIt != T[*tIt].bag.end(); sIt++)
            bag.push_back((int)*sIt);
        V_T.push_back(bag);
    }
    
    boost::graph_traits<TD_tree_dec_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(T); eIt != eEnd; eIt++){
        std::map<boost::graph_traits<TD_tree_dec_t>::vertex_descriptor, unsigned int>::iterator v, w;
        v = vertex_map.find(boost::source(*eIt, T));
        w = vertex_map.find(boost::target(*eIt, T));
        E_T.push_back(v->second);
        E_T.push_back(w->second);
    }
}


/* PREPROCESSING */

int sage_PP_MD(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::PP_MD(G, T, lb);

    make_sage_decomp(T, V_T, E_T);

    return lb;
}


int sage_PP_FI_TM(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::PP_FI_TM(G, T, lb);

    make_sage_decomp(T, V_T, E_T);

    return lb;
}


/* LOWER BOUNDS */


int sage_deltaC_least_c(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::deltaC_least_c(G);
}

int sage_LBN_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::LBN_deltaC(G);
}

int sage_LBNC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::LBNC_deltaC(G);
}

int sage_LBP_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::LBP_deltaC(G);
}

int sage_LBPC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::LBPC_deltaC(G);
}


/* EXACT TREE DECOMPOSITIONS */

int sage_exact_decomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::exact_decomposition_cutset(G, T, lb);

    make_sage_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

int sage_exact_decomposition_chordal(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::exact_decomposition_chordal(G, T);

    make_sage_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

/* APPOXIMATIVE TREE DECOMPOSITIONS */


int sage_seperator_algorithm(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::seperator_algorithm(G, T);

    make_sage_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

int sage_ordering_to_treedec(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, std::vector<unsigned int> &elim_ordering){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;
    treedec::ordering_to_treedec(G, elim_ordering, T);

    make_sage_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

void sage_treedec_to_ordering(std::vector<std::vector<int> > &V, std::vector<unsigned int> &E, std::vector<unsigned int> &elim_ordering){
    TD_tree_dec_t T;
    make_tdlib_decomp(T, V, E);

    treedec::treedec_to_ordering(T, elim_ordering);
}


/* MISC */


int sage_is_valid_decomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;
    make_tdlib_decomp(T, V_T, E_T);

    return treedec::is_valid_treedecomposition(G, T);
}

