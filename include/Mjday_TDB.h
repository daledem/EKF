//
// Created by David on 28/04/2024.
//

#ifndef MJD_TDB_H
#define MJD_TDB_H

#include <cmath>

//--------------------------------------------------------------------------
//
// Mjday_TDB: Computes the Modified Julian Date for barycentric dynamical
//            time
//
//  Inputs:
//    Mjd_TT      - Modified julian date (TT)
//
//  Output:
//    Mjd_TDB     - Modified julian date (TDB)
//
// Reference:
// Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill;
// New York; 3rd edition(2007).
//
// Last modified:   2015/08/12   M. Mahooti
//
//--------------------------------------------------------------------------

double Mjday_TDB(double Mjd_TT);

#endif //MJD_TDB_H
