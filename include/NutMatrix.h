#ifndef NUTMATRIX_H
#define NUTMATRIX_H

#include "./Matrix.h"
#include "./MeanObliquity.h"
#include "./NutAngles.h"
#include "./R_x.h"
#include "./R_y.h"
#include "./R_z.h"

//------------------------------------------------------------------------------
// NutMatrix(double Mjd_TT)
//------------------------------------------------------------------------------
/**
*   Transformation from mean to true equator and equinox
*
* @param <Mjd_TT> Modified Julian Date (Terrestrial Time)
*
* @return Nutation matrix
*
*/
//------------------------------------------------------------------------------

Matrix NutMatrix(double Mjd_TT);

#endif //NUTMATRIX_H
