//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/10
//
//------------------------------------------------------------------------------
#include "../include/NutMatrix.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix NutMatrix(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 *   Transformation from mean to true equator and equinox
 *
 * @param <Mjd_TT> Modified Julian Date (Terrestrial Time)
 *
 * @return Nutation matrix
 *
 */
//------------------------------------------------------------------------------
Matrix NutMatrix(double Mjd_TT) {
    double eps, dpsi, deps;
    // Mean obliquity of the ecliptic
    eps = MeanObliquity (Mjd_TT);

    // Nutation in longitude and obliquity
    NutAngles (dpsi, deps,Mjd_TT);

    // Transformation from mean to true equator and equinox
    return R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);

}

