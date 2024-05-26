//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/11
//
//------------------------------------------------------------------------------
#include "../include/LTC.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix LTC(double lon, double lat)
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
Matrix LTC(double lon, double lat) {
    double Aux;
    Matrix M = R_y(-1.0*lat)*R_z(lon);

    for(int j = 1; j <= 3; j++) {
        Aux=M(1,j);
        M(1,j)=M(2,j);
        M(2,j)=M(3,j);
        M(3,j)= Aux;
    }

    return M;
}
