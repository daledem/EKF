//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/08
//
//------------------------------------------------------------------------------
#include "../include/gast.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// double gast(double Mjd_UT1)
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
double gast(double Mjd_UT1) {
    return fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*Const::pi );
}
