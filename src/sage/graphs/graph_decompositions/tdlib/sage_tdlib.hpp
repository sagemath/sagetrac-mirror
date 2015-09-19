/* TdLib interface for sage */


/* PREPROCESSING */

int sage_preprocessing(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &bags, int lb);
int sage_preprocessing_MD(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb);
int sage_preprocessing_FI_TM(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb);

/* LOWER BOUNDS */

int sage_deltaC_min_d(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);
int sage_deltaC_max_d(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);
int sage_deltaC_least_c(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);

int sage_LBN_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);
int sage_LBNC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);
int sage_LBP_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);
int sage_LBPC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);


/* EXACT TREE DECOMPOSITIONS */

int sage_CR_greedy_decomp(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb);
int sage_CR_dynamic_decomp(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb);


/* APPROXIMATIVE TREE DECOMPOSITIONS */

int sage_seperator_algorithm(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);


/* POSTPROCESSING */

int sage_MSVS(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);
void sage_minimalChordal(std::vector<unsigned int> &V, std::vector<unsigned int> &E, std::vector<unsigned int> &old_elimination_ordering, std::vector<unsigned int> &new_elimination_ordering);


/* MISC */

int sage_is_valid_decomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);
int sage_ordering_to_treedec(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, std::vector<unsigned int> &elim_ordering);
void sage_treedec_to_ordering(std::vector<std::vector<int> > &V, std::vector<unsigned int> &E, std::vector<unsigned int> &elim_ordering);
int sage_get_width(std::vector<std::vector<int> > &V_T);

