#ifndef POLE_H
#define POLE_H

#include "./Matrix.h"
#include "./R_x.h"
#include "./R_y.h"

//------------------------------------------------------------------------------
// NutMatrix(double Mjd_TT)
//------------------------------------------------------------------------------
/**
*   Transformation from pseudo Earth-fixed to Earth-fixed coordinates
*      for a given date
*
* @param <xp> Pole coordinte
* @param <yp> Pole coordinte
*
* @return Pole matrix
*
*/
//------------------------------------------------------------------------------

Matrix PoleMatrix(double xp,double yp);

#endif //POLE_H
