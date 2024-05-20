#ifndef POSITION_H
#define POSITION_H

#include "./Matrix.h"
#include "./SAT_Const.h"
#include <cmath>

//------------------------------------------------------------------------------
// Position(double lon, double lat, double h)
//------------------------------------------------------------------------------
/**
*   Position vector (r [m]) from geodetic coordinates (Longitude [rad],
*       latitude [rad], altitude [m])
*
* @param <lon> longitude
* @param <lat> latitude
* @param <h> altitude
*
* @return Position vector ([m]) from geodetic coordinates
*
*/
//------------------------------------------------------------------------------

Matrix Position(double lon, double lat, double h);

#endif //POSITION_H
