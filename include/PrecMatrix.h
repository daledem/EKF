#ifndef PRECMATRIX_H
#define PRECMATRIX_H

#include "./Matrix.h"
#include "./R_z.h"
#include "./R_y.h"
#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// PrecMatrix(double Mjd_1,double Mjd_2)
//------------------------------------------------------------------------------
/**
*   Precession transformation of equatorial coordinates
*
* @param <Mjd_1> Epoch given (Modified Julian Date TT)
* @param <MjD_2> Epoch to precess to (Modified Julian Date TT)
*
* @return Precession transformation matrix
*
*/
//------------------------------------------------------------------------------

Matrix PrecMatrix(double Mjd_1,double Mjd_2);

#endif //PRECMATRIX_H
