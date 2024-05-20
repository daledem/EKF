#ifndef GHAMATRIX_H
#define GHAMATRIX_H

#include "./Matrix.h"
#include "./R_z.h"
#include "./gast.h"

//------------------------------------------------------------------------------
// GHAMatrix(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
*   Transformation from true equator and equinox to Earth equator and
*   Greenwich meridian system
*
* @param <Mjd_UT1> Modified Julian Date UT1
*
* @return Greenwich Hour Angle matrix
*
*/
//------------------------------------------------------------------------------

Matrix GHAMatrix(double Mjd_UT1);

#endif //GHAMATRIX_H
