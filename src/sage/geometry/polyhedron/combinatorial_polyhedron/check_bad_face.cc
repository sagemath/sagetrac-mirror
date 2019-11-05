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


inline uint64_t bit_lookup(uint64_t F, size_t n){
    // Return True if F contains bit corresponding to n.
    uint64_t tmp = 1;
    tmp = tmp << (64 - n - 1);
    return tmp & F;
}

inline vector<MachineInteger> make_ineq(size_t *PolyIneq, size_t mult){
    // make Hrep inequality
    vector<MachineInteger> ineq(mult);
    ineq[PolyIneq[0]] += 1;
    ineq[PolyIneq[1]] += 1;
    ineq[PolyIneq[2]] = -1;
    if (PolyIneq[0] + PolyIneq[1] + 2 > mult - 1)
        ineq[mult -1] += 1;
    return ineq;
}

int check_bad_face(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, size_t *Hrep, size_t n_Hrep, size_t e){
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
            FaceEq.append(make_ineq(PolyIneq[j],mult));
            if (counter < n_Hrep-1){
                counter++;
            }
        }
        else
            FaceInEq.append(make_ineq(PolyIneq[j],mult));
    }

    Matrix<MachineInteger> Frob(0,mult);

    for(size_t f0=0;f0<dim;f0++){  //f will run from 1
        if (bit_lookup(FPos, f0))
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
        WilfNeg[dim]=f-mult-e*(f-mult)-e;
        Frob.append(WilfNeg);

        Cone<MachineInteger> WilfPolyhedron(Type::inhom_inequalities, Frob,
                                            Type::inhom_equations, FaceEq, Type::strict_signs, StrictSigns);
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
