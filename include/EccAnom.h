#ifndef ECCANOM_H
#define ECCANOM_H

#include <cmath>
#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// EccAnom(double M, double e)
//------------------------------------------------------------------------------
/**
*   Computes the eccentric anomaly for elliptic orbits
*
* @param <M> Mean anomaly in [rad]
* @param <e> Eccentricity of the orbit [0,1]
*
* @return Eccentric anomaly in [rad]
*/
//------------------------------------------------------------------------------

double EccAnom(double M, double e);

#endif //ECCANOM_H
