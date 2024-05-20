//
// Created by David on 28/04/2024.
//

#ifndef MJD_TDB_H
#define MJD_TDB_H

#include <cmath>

//------------------------------------------------------------------------------
// Mjday_TDB(double Mjd_TT)
//------------------------------------------------------------------------------
/**
*   Computes the Modified Julian Date for barycentric dynamical
*       time
*
* @param <Mjd_TT> Modified julian date (TT)
*
* @return Modified julian date (TDB)
*
*/
//------------------------------------------------------------------------------

double Mjday_TDB(double Mjd_TT);

#endif //MJD_TDB_H
