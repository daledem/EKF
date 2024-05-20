#ifndef MEANOBLIQUITY_H
#define MEANOBLIQUITY_H

#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// MeanObliquity(double Mjd_TT)
//------------------------------------------------------------------------------
/**
*   Computes the mean obliquity of the ecliptic
*
* @param <Mjd_TT> Geodetic East longitude [rad]
*
* @return Mean obliquity of the ecliptic [rad]
*
*/
//------------------------------------------------------------------------------

double MeanObliquity(double Mjd_TT);

#endif //MEANOBLIQUITY_H
