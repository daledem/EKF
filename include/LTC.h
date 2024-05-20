#ifndef LTC_H
#define LTC_H

#include "./Matrix.h"
#include "./R_y.h"
#include "./R_z.h"

//------------------------------------------------------------------------------
// LTC(double lon, double lat)
//------------------------------------------------------------------------------
/**
*   Transformation from Greenwich meridian system to
*      local tangent coordinates
*
* @param <lon> Geodetic East longitude [rad]
* @param <lat> Geodetic latitude [rad]
*
* @return Rotation matrix from the Earth equator and Greenwich meridian
*             to the local tangent (East-North-Zenith) coordinate system
*
*/
//------------------------------------------------------------------------------

Matrix LTC(double lon, double lat);

#endif //LTC_H
