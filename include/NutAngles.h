#ifndef NUTANGLES_H
#define NUTANGLES_H

#include "./Matrix.h"
#include "./SAT_Const.h"
#include <cmath>

//------------------------------------------------------------------------------
// NutAngles(double& dpsi,double& deps,double Mjd_TT)
//------------------------------------------------------------------------------
/**
*   Nutation in longitude and obliquity
*
* @param <dpsi> Nutation Angles
* @param <deps> Nutation Angles
* @param <Mjd_TT> Modified Julian Date (Terrestrial Time)
*
*/
//------------------------------------------------------------------------------

void NutAngles(double& dpsi,double& deps,double Mjd_TT);

#endif //NUTANGLES_H
