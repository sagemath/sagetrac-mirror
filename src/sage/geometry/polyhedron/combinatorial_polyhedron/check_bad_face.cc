/*
Much of the code is copied from https://github.com/Normaliz/Normaliz/blob/wilf/source/libnormaliz/cone.cpp

We refer to "Wilf's conjecture in fixed multiplicity" by
Bruns, Garcia-Sanchez, O'Neill and Wilburne as [BGOW2019].

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
    // Return True if g*F contains bit corresponding to n.
    // g ist the group action, of which
    // ``group_action_inv` is the invers.
    uint64_t tmp = 1;

    // Apply ``group_action_inv`` to ``n``.
    size_t new_n = (n*group_action_inv + group_action_inv - 1) % mult;

    //size_t new_n = n;
    tmp = tmp << (64 - new_n - 1);
    return tmp & F;
}

inline void addvector(vector<MachineInteger> *a, vector<MachineInteger> *b, MachineInteger times){
    // Set a += times*b
    size_t len = a[0].size();
    for(size_t i=0; i< len; i++){
        a[0][i] += times*b[0][i];
    }
}

inline vector<MachineInteger> make_eq(size_t *PolyIneq, size_t mult, size_t group_action){
    // Return the defining equality induced by
    // x_{g*i} and x_{g*j} on the LHS and
    // x_{g*(i+j)} on the RHS,
    // where g denotes the group action in (ZZ/mult*ZZ)^*
    // and i,j,i+j the elements in PolyIneq[0..2] resp.
    //
    // See Definition 3.1 in [BGOW2019].
    //
    // The equality is given by
    // x_{g*i} + x_{g*j} = x_{g*(i+j)}, if (g*i)+(g*j) < m
    // and
    // x_{g*i} + x_{g*j} + 1 = x_{g*(i+j)}, if (g*i)+(g*j) > m
    //
    // The output is a vector with the indices of the
    // basis elements x_1,...,x_{mult},1.
    vector<MachineInteger> eq(mult);

    // Applying the group action.
    size_t first_left   = (PolyIneq[0]*group_action + group_action - 1) % mult;
    size_t second_left  = (PolyIneq[1]*group_action + group_action - 1) % mult;
    size_t right        = (PolyIneq[2]*group_action + group_action - 1) % mult;

    eq[first_left] += 1;
    eq[second_left] += 1;
    eq[right] = -1;

    if (first_left + second_left + 2 > mult - 1)
        eq[mult -1] += 1;
    return eq;
}

inline vector<MachineInteger> make_ineq(size_t *PolyIneq, size_t mult, size_t group_action){
    // Return the defining inequality induced by
    // x_{g*i} and x_{g*j} on the LHS and
    // x_{g*(i+j)} on the RHS
    // for the interior of the face,
    // where g denotes the group action in (ZZ/mult*ZZ)^*
    // and i,j,i+j the elements in PolyIneq[0..2] resp.
    //
    // See Definition 3.1 in [BGOW2019].
    //
    // The inequality is given by
    // x_{g*i} + x_{g*j} > x_{g*(i+j)}, if (g*i)+(g*j) < m
    // and
    // x_{g*i} + x_{g*j} + 1 > x_{g*(i+j)}, if (g*i)+(g*j) > m
    //
    // or equivalently (only integer interior points)
    //
    // x_{g*i} + x_{g*j} >= x_{g*(i+j)} + 1, if (g*i)+(g*j) < m
    // and
    // x_{g*i} + x_{g*j} + 1 >= x_{g*(i+j)} + 1, if (g*i)+(g*j) > m
    //
    // The output is a vector with the indices of the
    // basis elements x_1,...,x_{mult},1.
    //
    // NOTE: This is the same as the corresponding
    // inequality, but a -1 added to the LHS.
    vector<MachineInteger> ineq = make_eq(PolyIneq, mult, group_action);
    ineq[mult -1] -= 1;
    return ineq;
}

inline vector<MachineInteger> make_eq(size_t *PolyIneq, size_t len, Matrix<MachineInteger> *Basis, \
                                        size_t mult, size_t group_action){
    // Return the defining equality with respect to a
    // new basis,
    // where the first element in the new Basis
    // corresponds to x_1, the second to x_2 and so on.
    vector<MachineInteger> eq(len);

    // Applying the group action.
    size_t first_left   = (PolyIneq[0]*group_action + group_action - 1) % mult;
    size_t second_left  = (PolyIneq[1]*group_action + group_action - 1) % mult;
    size_t right        = (PolyIneq[2]*group_action + group_action - 1) % mult;
    addvector(&eq, &Basis[0][first_left], 1);
    addvector(&eq, &Basis[0][second_left], 1);
    addvector(&eq, &Basis[0][right], -1);

    if (first_left + second_left + 2 > mult - 1)
        eq[len-1] += 1;
    return eq;
}

inline vector<MachineInteger> make_ineq(size_t *PolyIneq, size_t len, Matrix<MachineInteger> *Basis, \
                                        size_t mult, size_t group_action){
    // Return the defining inequality with respect to a
    // new basis,
    // where the first element in the new Basis
    // corresponds to x_1, the second to x_2 and so on.
    vector<MachineInteger> ineq = make_eq(PolyIneq, len, Basis, mult, group_action);
    ineq[len -1] -= 1;
    return ineq;
}

inline MachineInteger maxpos(vector<MachineInteger> *a, size_t *index){
    // Return the largest positive entry of ``a`` and set
    // index to the corresponding entry.
    // Return -1, if there is no positive entry
    // and at least one negative entry.
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
    // Find a positive value of ``Game``,
    // where the corresponding entry is not negative in
    // ``a``.
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
    // See if adding b to a will do any good.
    //
    // We try to add vectors to ``a`` such that ``a``
    // becomes non-positive (and at least one negative).
    //
    // Every vector ``b`` will get a score, according to
    // who good we think it is.
    //
    // We will get points for
    // - subtracting from positive
    // - negative for adding to positive
    // (normalize the first two by the total number of
    // positives)
    // and negative points for
    // - adding to negative
    // (normalize by the total of negatives).
    //
    // However, we will not make something negative
    // positive again ever.
    // It seems that the negative entries in the negation
    // of (4.2) see [BGOW2019] are really hard to get
    // negative. As we start with the negation of (4.2)
    // it seems pointless to make negative entries
    // positive again.
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

inline int check_bad_face(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, uint8_t *Hrep, size_t n_Hrep, size_t e, size_t group_action, size_t group_action_inv){
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

    if (verbose){
        cout << "Length of Hrep " << n_Hrep << endl;
    }

    Matrix<MachineInteger> Basis(0,mult);
    // Make some basis vectors.
    for(size_t j=0;j<mult-1;++j){
        vector<MachineInteger> bas(mult);
        bas[j] = 1;
        Basis.append(bas);
    }

    // Make basis substition based on equalities.
    size_t length_new_basis = mult;
    for (size_t j=0;j<n_Hrep;j++){
        vector<MachineInteger> equality = make_eq(PolyIneq[Hrep[j]],mult,&Basis, mult, group_action);
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
    Matrix<MachineInteger> Signs(0,length_new_basis);
    // Make some basis vectors.
    for(size_t j=0;j<mult-1;++j){
        vector<MachineInteger> bas(length_new_basis);
        vector<MachineInteger> signs(length_new_basis);
        counter = 0;
        for(size_t i=0;i<mult-1;i++){
            if (Basis[i][i] == 1){
               bas[counter] = Basis[j][i];
               counter++;
            }
        }
        bas[length_new_basis -1] = Basis[j][mult-1];
        addvector(&signs, &bas, 1);
        signs[length_new_basis-1] += -1;

        Basis2.append(bas);
        Signs.append(signs);
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
            FaceInEq2.append(make_ineq(PolyIneq[j],length_new_basis,&Basis2, mult, group_action));
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
            addvector(&frob, &Basis2[i], -1);
            if(i>f0)
                frob[length_new_basis -1] -= 1;
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
        WilfNeg[length_new_basis - 1] += f-mult-e*(f-mult)-e-1;
        Frob.append(WilfNeg);

        // Trying to play a game instead of creating the
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
        // We delete all the vectors that are strictly
        // elementwise larger than others.
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
        for(size_t j=0; j< 200; j++){
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

        Frob2.append(Signs);

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
            continue;
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

inline int check_bad_face(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, size_t *Hrep, size_t n_Hrep, size_t e, size_t group_action, size_t group_action_inv){
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
    for(size_t j=0;j<mult-1;++j){
        vector<MachineInteger> bas(mult);
        bas[j] = 1;
        Basis.append(bas);
    }

    // Make basis substition based on equalities.
    size_t length_new_basis = mult;
    for (size_t j=0;j<n_Hrep;j++){
        vector<MachineInteger> equality = make_eq(PolyIneq[Hrep[j]],mult,&Basis, mult, group_action);
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
    Matrix<MachineInteger> Signs(0,length_new_basis);
    // Make some basis vectors.
    for(size_t j=0;j<mult-1;++j){
        vector<MachineInteger> bas(length_new_basis);
        vector<MachineInteger> signs(length_new_basis);
        counter = 0;
        for(size_t i=0;i<mult-1;i++){
            if (Basis[i][i] == 1){
               bas[counter] = Basis[j][i];
               counter++;
            }
        }
        bas[length_new_basis -1] = Basis[j][mult-1];
        addvector(&signs, &bas, 1);
        signs[length_new_basis-1] += -1;

        Basis2.append(bas);
        Signs.append(signs);
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
            FaceInEq2.append(make_ineq(PolyIneq[j],length_new_basis,&Basis2, mult, group_action));
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
            addvector(&frob, &Basis2[i], -1);
            if(i>f0)
                frob[length_new_basis -1] -= 1;
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
        WilfNeg[length_new_basis - 1] += f-mult-e*(f-mult)-e-1;
        Frob.append(WilfNeg);

        // Trying to play a game instead of creating the
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
        // We delete all the vectors that are strictly
        // elementwise larger than others.
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
        for(size_t j=0; j< 200; j++){
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

        Frob2.append(Signs);

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
            continue;
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

inline int check_bad_face_original(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, uint8_t *Hrep, size_t n_Hrep, size_t e, size_t group_action, size_t group_action_inv){
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

    vector<MachineInteger> all_one(dim,1);
    Matrix<MachineInteger> StrictSigns(1,dim);
    StrictSigns[0]=all_one;

    Matrix<MachineInteger> FaceEq(0,mult);
    Matrix<MachineInteger> FaceInEq(0,mult);

    size_t counter = 0;
    for(size_t j=0;j<n_coatoms;++j){

        if(Hrep[counter] == j){
            FaceEq.append(make_eq(PolyIneq[j],mult, group_action));
            if (counter < n_Hrep-1){
                counter++;
            }
        }
        else {
            FaceInEq.append(make_ineq(PolyIneq[j],mult, group_action));
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
                                            Type::inhom_equations, FaceEq,
                                            Type::strict_signs, StrictSigns);
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

inline int check_bad_face_original(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, size_t *Hrep, size_t n_Hrep, size_t e, size_t group_action, size_t group_action_inv){
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

    vector<MachineInteger> all_one(dim,1);
    Matrix<MachineInteger> StrictSigns(1,dim);
    StrictSigns[0]=all_one;

    Matrix<MachineInteger> FaceEq(0,mult);
    Matrix<MachineInteger> FaceInEq(0,mult);

    size_t counter = 0;
    for(size_t j=0;j<n_coatoms;++j){

        if(Hrep[counter] == j){
            FaceEq.append(make_eq(PolyIneq[j],mult, group_action));
            if (counter < n_Hrep-1){
                counter++;
            }
        }
        else {
            FaceInEq.append(make_ineq(PolyIneq[j],mult, group_action));
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
                                            Type::inhom_equations, FaceEq,
                                            Type::strict_signs, StrictSigns);
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
    // We are given a representative F of a bad orbit of the
    // KunzCone.
    //
    // We check, wether for every face in the orbit, there
    // is a counterexample to the Wilf conjecture (or
    // not).
    //
    // Return 0 if all such polyhedra are empty.
    // Return 1 if any such polyhedron is non-emtpy.
    // Return 2 if any such polyhedron contains an integer
    // point (yielding a counterexample to the Wilf
    // conjecture.
    int value = 0;
    for(size_t i=1; i < m; i++){
        size_t j=1;
        while (j < m){
            if ((i*j) % m == 1){
                break;
            }
            j++;
        }
        if ((i*j) % m != 1){
            continue;
        }
        int new_value = check_bad_face(PolyIneq, n_coatoms, m, LHS, Hrep, n_Hrep, e, i, j);
        if (new_value > value){
            value = new_value;
        }
        /*
        for(size_t j=1; j < m; j++){
            if ((i*j) % m == 1){
                // i and j are inverse in (ZZ/mZZ)*.
                // Hence, i*F is a face in the orbit.
                // Check, wether the conjecture holds for
                // i*F.
                int new_value = check_bad_face(PolyIneq, n_coatoms, m, LHS, Hrep, n_Hrep, e, i, j);
                if (new_value > value){
                    value = new_value;
                }
            }
        }
        */
    }
    return value;
}

int check_bad_faces(size_t **PolyIneq, size_t n_coatoms, size_t m, \
                    uint64_t *LHS, size_t n_bad_faces, uint8_t **bad_faces){
    // We are given a representative F of a bad orbit of the
    // KunzCone.
    //
    // We check, wether for every face in the orbit, there
    // is a counterexample to the Wilf conjecture (or
    // not).
    //
    // Return 0 if all such polyhedra are empty.
    // Return 1 if any such polyhedron is non-emtpy.
    // Return 2 if any such polyhedron contains an integer
    // point (yielding a counterexample to the Wilf
    // conjecture.
    int value = 0;
    for(size_t k=0; k<n_bad_faces; k++){
        for(size_t i=1; i < m; i++){
            size_t j=1;
            while (j < m){
                if ((i*j) % m == 1){
                    break;
                }
                j++;
            }
            if ((i*j) % m != 1){
                continue;
            }
            if (n_bad_faces == 0){
                return 0;
            }
            size_t e = (size_t) bad_faces[k][1];
            size_t n_Hrep = (size_t) bad_faces[k][0];
            int new_value = check_bad_face(PolyIneq, n_coatoms, m, LHS[k], bad_faces[k] + 2, n_Hrep, e, i, j);
            if (new_value > value){
                value = new_value;
            }
        }
        /*
        for(size_t j=1; j < m; j++){
            if ((i*j) % m == 1){
                // i and j are inverse in (ZZ/mZZ)*.
                // Hence, i*F is a face in the orbit.
                // Check, wether the conjecture holds for
                // i*F.
                int new_value = check_bad_face(PolyIneq, n_coatoms, m, LHS, Hrep, n_Hrep, e, i, j);
                if (new_value > value){
                    value = new_value;
                }
            }
        }
        */
    }
    return value;
}
