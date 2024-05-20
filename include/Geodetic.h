#ifndef GEODETIC_H
#define GEODETIC_H

#include <cmath>
#include "./Matrix.h"
#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// Geodetic(double& lon, double& lat, double& h, const Matrix& r)
//------------------------------------------------------------------------------
/**
*   geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
*   from given position vector (r [m])
*
* @param <lon> Longitude [rad]
* @param <lat> latitude [rad]
* @param <h> altitude [m]
* @param <r> position vector [m]
*
* @return Geodetic coordinates
*
*/
//------------------------------------------------------------------------------

void Geodetic(double& lon, double& lat, double& h, const Matrix& r);

#endif //GEODETIC_H
