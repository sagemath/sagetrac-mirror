/*
This code is mostly copied from https://github.com/Normaliz/Normaliz/blob/wilf/source/libnormaliz/cone.cpp

#*****************************************************************************
#
#       Copyright (C) 2007-2014 Winfried Bruns, Bogdan Ichim, Christof Soeger
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
*/

#include <stdlib.h>
#include <list>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <cstdint>
#include <cstdio>
using namespace std;
#include <libnormaliz/cone.h>
#include <libnormaliz/vector_operations.h>
#include <libnormaliz/map_operations.h>
#include <libnormaliz/convert.h>
#include <libnormaliz/my_omp.h>
using namespace libnormaliz;


inline uint64_t bit_lookup(uint64_t F, size_t n, size_t mult, size_t group_action_inv){
    // Return True if F contains bit corresponding to n.
    uint64_t tmp = 1;
    size_t new_n = (n*group_action_inv + group_action_inv - 1) % mult;
    //size_t new_n = n;
    tmp = tmp << (64 - new_n - 1);
    return tmp & F;
}

inline vector<MachineInteger> make_ineq(size_t *PolyIneq, size_t mult, int is_equality, size_t group_action){
    // make Hrep inequality
    vector<MachineInteger> ineq(mult);
    size_t first_left   = (PolyIneq[0]*group_action + group_action - 1) % mult;
    size_t second_left  = (PolyIneq[1]*group_action + group_action - 1) % mult;
    size_t right        = (PolyIneq[2]*group_action + group_action - 1) % mult;
    ineq[first_left] += 1;
    ineq[second_left] += 1;
    ineq[right] = -1;
    if (first_left + second_left + 2 > mult - 1)
        ineq[mult -1] += 1;
    if (!is_equality)
        ineq[mult -1] -= 1;
    return ineq;
}

inline void addvector(vector<MachineInteger> *a, vector<MachineInteger> *b, MachineInteger times){
    // Set a += times*b
    size_t len = a[0].size();
    for(size_t i=0; i< len; i++){
        a[0][i] += times*b[0][i];
    }
}

inline vector<MachineInteger> make_ineq(size_t *PolyIneq, size_t len, int is_equality, Matrix<MachineInteger> *Basis, size_t mult, size_t group_action){
    // make Hrep inequality
    // with Basis specified on input.
    vector<MachineInteger> ineq(len);
    size_t first_left   = (PolyIneq[0]*group_action + group_action - 1) % mult;
    size_t second_left  = (PolyIneq[1]*group_action + group_action - 1) % mult;
    size_t right        = (PolyIneq[2]*group_action + group_action - 1) % mult;
    addvector(&ineq, &Basis[0][first_left], 1);
    addvector(&ineq, &Basis[0][second_left], 1);
    addvector(&ineq, &Basis[0][right], -1);
    if (first_left + second_left + 2 > mult - 1)
        addvector(&ineq, &Basis[0][mult-1], 1);
    if (!is_equality)
        addvector(&ineq, &Basis[0][mult-1], -1);
    return ineq;
}

inline MachineInteger maxpos(vector<MachineInteger> *a, size_t *index){
    // Return the largest positive entry of a and set
    // index to the corresponding entry.
    // Return -1, if there is no positive entry.
    // (In this case we have already won our game.)
    MachineInteger output = 0;
    for(size_t i=0; i < a[0].size(); i++){
        if(a[0][i] > output){
            output = a[0][i];
            index[0] = i;
        }
    }
    if (output == 0){
        for(size_t i=0; i < a[0].size(); i++){
            if(a[0][i] < 0){
                return -1;
            }
        }
    }
    return output;
}

inline void find_unchanged_pos(vector<MachineInteger> *Game, vector<MachineInteger> *a, size_t *index){
    for(size_t i=0; i < a[0].size(); i++){
        if ((Game[0][i] > 0) && (a[0][i] >= 0)){
            index[0] = i;
            return;
        }
    }
}

inline int is_smaller(vector<MachineInteger> *a, vector<MachineInteger> *b){
    // Return 1 if a[0][i] <= b[0][i] for all i.
    for(size_t i=0; i< a[0].size(); i++){
        if(a[0][i] > b[0][i]){
            return 0;
        }
    }
    return 1;
}

inline MachineInteger combine_score(vector<MachineInteger> *a, vector<MachineInteger> *b){
    //See if adding b to a will do any good.
    MachineInteger prevpos=0, prevneg=0, pos_to_pos=0, pos_to_neg=0, neg_to_pos=0, neg_to_neg=0, sum=0;
    for(size_t i=0; i<a[0].size();i++){
        if(a[0][i] > 0){
            prevpos += a[0][i];
            if(b[0][i] > 0){
                sum += b[0][i];
                pos_to_pos += b[0][i];
            } else {
                sum -= b[0][i];
                neg_to_pos -= b[0][i];
            }
        } else {
            prevneg -= a[0][i];
            if(b[0][i] > 0){
                if (b[0][i] > -a[0][i])
                    // Do not make something negative
                    // positiv again.
                    return -1000000;
                sum += b[0][i];
                pos_to_neg += b[0][i];
            } else {
                neg_to_neg -= b[0][i];
                sum -= b[0][i];
            }
        }
    }
    // Trying to make sense of it all.
    if (prevpos == 0){
        // Looks like we have won already.
        return -1000000;
    }
    if (prevneg == 0){
        // We have lost anyway already.
        return -1000000;
    }
    if (sum == 0){
        return -1000;
    }
    return ((neg_to_pos - pos_to_pos)*10000/prevpos - ((pos_to_neg)*10000/prevneg))/sum;
}

int check_bad_face(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, size_t *Hrep, size_t n_Hrep, size_t e, size_t group_action, size_t group_action_inv){
    // For a bad face in in the Kunz Cone we check
    // wether with the inequalites from (4.1) and negation
    // of (4.2)
    // (see paper "Wilf's conjecture in fixed multiplicity")
    // yield a nonempty polyhedron and wether it contains
    // an interior point.
    //
    // Return 0 if all such polyhedra are empty.
    // Return 1 if any such polyhedron is non-emtpy.
    // Return 2 if any such polyhedron contains an integer
    // point (yielding a counterexample to the Wilf
    // conjecture.

    uint64_t FPos = LHS;
    size_t dim = m-1;
    size_t mult = m;
    int output = 0;
    int verbose = 0;

    Matrix<MachineInteger> Basis(0,mult);
    // Make some basis vectors.
    for(size_t j=0;j<mult;++j){
        vector<MachineInteger> bas(mult);
        bas[j] = 1;
        Basis.append(bas);
    }

    // Make basis substition based on equalities.
    size_t length_new_basis = mult;
    for (size_t j=0;j<n_Hrep;j++){
        vector<MachineInteger> equality = make_ineq(PolyIneq[Hrep[j]],mult,1,&Basis, mult, group_action);
        for(size_t i=0;i<mult-1;i++){
            if (equality[i] == 1){
                length_new_basis -= 1;
                for(size_t j=0; j<mult-1; j++){
                    addvector(&Basis[j], &equality, -Basis[j][i]);
                }
                break;
            }
            else if (equality[i] == -1){
                length_new_basis -= 1;
                for(size_t j=0; j<mult-1; j++){
                    addvector(&Basis[j], &equality, Basis[j][i]);
                }
                break;
            }
        }
    }
    if (verbose){
        cout << "This is the new basis " << endl;
        Basis.pretty_print(cout);
    }

    size_t counter = 0;
    Matrix<MachineInteger> Basis2(0,length_new_basis);
    // Make some basis vectors.
    for(size_t j=0;j<mult;++j){
        vector<MachineInteger> bas(length_new_basis);
        counter = 0;
        for(size_t i=0;i<mult;i++){
            if (Basis[i][i] == 1){
               bas[counter] = Basis[j][i];
               counter++;
            }
        }
        Basis2.append(bas);
    }
    if (verbose){
        cout << "This is the transformed basis " << endl;
        Basis2.pretty_print(cout);
    }

    Matrix<MachineInteger> FaceInEq2(0,length_new_basis);
    counter = 0;
    for(size_t j=0;j<n_coatoms;++j){

        if(Hrep[counter] == j){
            if (counter < n_Hrep-1){
                counter++;
            }
        }
        else {
            FaceInEq2.append(make_ineq(PolyIneq[j],length_new_basis,0,&Basis2, mult, group_action));
        }
    }
    if (verbose){
        cout << "Those are the new inequalites " << endl;
        FaceInEq2.pretty_print(cout);
    }


    Matrix<MachineInteger> Frob(0,length_new_basis);

    for(size_t f0=0;f0<dim;f0++){  //f will run from 1
        if (bit_lookup(FPos, f0, mult, group_action_inv))
                continue;

        // now we are in business
        Frob.resize(0,true);

        // make Frobenius number inequalities
        for(size_t i=0;i<dim;++i){
            if(i==f0)
                continue;
            vector<MachineInteger> frob(length_new_basis);
            addvector(&frob, &Basis2[f0], 1);
            //frob[f0]=1;
            addvector(&frob, &Basis2[i], -1);
            //frob[i]=-1;
            if(i>f0)
                addvector(&frob, &Basis2[dim], -1);
            Frob.append(frob);
        }

        Frob.append(FaceInEq2); // combining inequalities, all collected in Frob
        // now the Wilf negation

        vector<MachineInteger> WilfNeg(length_new_basis);
        for(size_t i=0;i<dim;++i){
            if(i==f0)
                addvector(&WilfNeg, &Basis2[i], mult-e*mult+e);
                //WilfNeg[i]=mult-e*mult+e;
            else
                addvector(&WilfNeg, &Basis2[i], e);
                //WilfNeg[i]=e;
        }
        long f=f0+1;
        addvector(&WilfNeg, &Basis2[dim], f-mult-e*(f-mult)-e-1);
        Frob.append(WilfNeg);

        // Trying to play the game instead of creating the
        // cone.
        vector<MachineInteger> Game(length_new_basis);
        addvector(&Game, &WilfNeg, 1);
        if (verbose)
            cout << "Start game with " << Game << endl;

        size_t index = 0;
        size_t maxindex = 0;
        index = 0;
        MachineInteger maximum = 0;
        MachineInteger score;
        MachineInteger best_score;
        maximum = maxpos(&Game, &maxindex);
        if ((maximum == -1)){
            // We have won the game.
            if (verbose)
                cout << "Have won the game with " << Game << endl;
            continue;
        }

        // Creating a cleaned up version of Frob.
        // We delete all the vectors that are trivially
        Matrix<MachineInteger> Frob2(0,length_new_basis);

        size_t counter = 0;
        for(size_t i=0; i< dim + n_coatoms - n_Hrep-1; i++){
            int add = 1;
            for(size_t j=0; j < counter; j++){
                if (is_smaller(&Frob[i], &Frob2[j])){
                    add = 0;
                    Frob2[j] = Frob[i];
                    break;
                }
                else {
                    if (is_smaller(&Frob2[j], &Frob[i])){
                        add = 0;
                        break;
                    }
                }
            }
            if (add){
                Frob2.append(Frob[i]);
                counter++;
            }
        }
        Frob2.append(Frob[dim + n_coatoms - n_Hrep-1]);

        // How many rounds we try.
        for(size_t j=0; j< 1000; j++){
            best_score = -1000000;
            for(size_t i=0; i< counter; i++){
                // We are looking for a inequality, that
                // will lower the maximum entry of our Game
                // inequality.
                if (Frob2[i][maxindex] < 0){
                    score = combine_score(&Game, &Frob2[i]);
                    if (score > best_score){
                        best_score = score;
                        index = i;
                    }
                }
            }
            if(best_score == -1000000){
                // Can't find an inequality anymore that
                // will help.
                // It seems our heuristics do not work in
                // this case.
                break;
            }
            addvector(&Game, &Frob2[index], 1);
            if (verbose){
                cout << "add             " << Frob2[index] << endl;
                cout << "to obtain       " << Game << endl;
            }
            maximum = maxpos(&Game, &maxindex);
            if ((maximum == -1)){
                // We have won the game.
                if (verbose)
                    cout << "Have won the game with " << Game << endl;
                break;
            }
            // Maybe there is a positive entry that we
            // have not touched (or even raised).
            // In this case we force the
            // next run to lower that index.
            find_unchanged_pos(&Game, &Frob2[index], &maxindex);
        }
        if (maximum == -1){
            // Go to next f.
            continue;
        }
        if (verbose)
            cout << "Could not win game with " << Game << endl;

        Cone<MachineInteger> WilfPolyhedron(Type::inhom_inequalities, Frob2);
        WilfPolyhedron.setVerbose(false);
        WilfPolyhedron.compute(ConeProperty::AffineDim);
        if(WilfPolyhedron.getAffineDim()>=0){
            cout << "Wilf polyhedron not empty" << endl;
            Frob.pretty_print(cout);
            cout << "----------------" << endl;
            output = 1;

        }
        else{
            if (verbose){
                cout << "\nHad to use the normaliz." << endl;
                Frob2.pretty_print(cout);
                cout << "----------------" << endl;
                Basis.pretty_print(cout);
                cout << "----------------" << endl;
                output = 0;
            }
        }

        if(WilfPolyhedron.getNrModuleGenerators()>0){
            cout << WilfPolyhedron.getNrModuleGenerators() << " counterexamples" << endl;
            WilfPolyhedron.getModuleGeneratorsMatrix().pretty_print(cout);
            cout << "f " << f << " e " << e << " mult " << mult << endl;
            return 2;
        }
    }
    return output;
}

int check_bad_face_original(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, size_t *Hrep, size_t n_Hrep, size_t e, size_t group_action, size_t group_action_inv){
    // For a bad face in in the Kunz Cone we check
    // wether with the inequalites from (4.1) and negation
    // of (4.2)
    // (see paper "Wilf's conjecture in fixed multiplicity")
    // yield a nonempty polyhedron and wether it contains
    // an interior point.
    //
    // Return 0 if all such polyhedra are empty.
    // Return 1 if any such polyhedron is non-emtpy.
    // Return 2 if any such polyhedron contains an integer
    // point (yielding a counterexample to the Wilf
    // conjecture.

    uint64_t FPos = LHS;
    size_t dim = m-1;
    size_t mult = m;
    int output = 0;

    Matrix<MachineInteger> FaceEq(0,mult);
    Matrix<MachineInteger> FaceInEq(0,mult);

    size_t counter = 0;
    for(size_t j=0;j<n_coatoms;++j){

        if(Hrep[counter] == j){
            FaceEq.append(make_ineq(PolyIneq[j],mult,1, group_action));
            if (counter < n_Hrep-1){
                counter++;
            }
        }
        else {
            FaceInEq.append(make_ineq(PolyIneq[j],mult,0, group_action));
        }
    }

    Matrix<MachineInteger> Frob(0,mult);

    for(size_t f0=0;f0<dim;f0++){  //f will run from 1
        if (bit_lookup(FPos, f0, mult, group_action_inv))
                continue;

        // now we are in business
        Frob.resize(0,true);

        // make Frobenius number inequalities
        for(size_t i=0;i<dim;++i){
            if(i==f0)
                continue;
            vector<MachineInteger> frob(mult);
            frob[f0]=1;
            frob[i]=-1;
            if(i>f0)
                frob[dim]=-1;
            Frob.append(frob);
        }

        Frob.append(FaceInEq); // combining inequalities, all collected in Frob
        // now the Wilf negation

        vector<MachineInteger> WilfNeg(mult);
        for(size_t i=0;i<dim;++i){
            if(i==f0)
                WilfNeg[i]=mult-e*mult+e;
            else
                WilfNeg[i]=e;
        }
        long f=f0+1;
        WilfNeg[dim]=f-mult-e*(f-mult)-e-1;
        Frob.append(WilfNeg);

        Cone<MachineInteger> WilfPolyhedron(Type::inhom_inequalities, Frob,
                                            Type::inhom_equations, FaceEq);
        WilfPolyhedron.setVerbose(false);
        WilfPolyhedron.compute(ConeProperty::ModuleGenerators);
        if(WilfPolyhedron.getAffineDim()>=0){
            cout << "Wilf polyhedron not empty" << endl;
            Frob.pretty_print(cout);
            cout << "----------------" << endl;
            FaceEq.pretty_print(cout);
            cout << "----------------" << endl;
            output = 1;

        }

        if(WilfPolyhedron.getNrModuleGenerators()>0){
            cout << WilfPolyhedron.getNrModuleGenerators() << " counterexamples" << endl;
            WilfPolyhedron.getModuleGeneratorsMatrix().pretty_print(cout);
            cout << "f " << f << " e " << e << " mult " << mult << endl;
            return 2;
        }
    }
    return output;
}

int check_bad_face(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, size_t *Hrep, size_t n_Hrep, size_t e){
    for(size_t i=1; i < m; i++){
        for(size_t j=1; j < m; j++){
            if ((i*j) % m == 1){
                int value = check_bad_face(PolyIneq, n_coatoms, m, LHS, Hrep, n_Hrep, e, i, j);
                if (value){
                    return value;
                }
            }
        }
    }
    return 0;
}

