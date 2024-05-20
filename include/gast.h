#ifndef GAST_H
#define GAST_H

#include <cmath>
#include "./gmst.h"
#include "./EqnEquinox.h"
#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// gast(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
*   Greenwich Apparent Sidereal Time
*
* @param <Mjd_UT1> Modified Julian Date UT1
*
* @return GAST in [rad]
*
*/
//------------------------------------------------------------------------------

double gast(double Mjd_UT1);

#endif //GAST_H
