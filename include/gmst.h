#ifndef GMST_H
#define GMST_H

#include <cmath>
#include "./Frac.h"
#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// gmst(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
*   Greenwich Mean Sidereal Time
*
* @param <Mjd_UT1> Modified Julian Date UT1
*
* @return GMST in [rad]
*
*/
//------------------------------------------------------------------------------

double gmst(double Mjd_UT1);

#endif //GMST_H
