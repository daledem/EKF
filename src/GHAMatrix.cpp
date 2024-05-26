//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/09
//
//------------------------------------------------------------------------------
#include "../include/GHAMatrix.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix GHAMatrix(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 *   Transformation from true equator and equinox to Earth equator and
 *   Greenwich meridian system
 *
 * @param <Mjd_UT1> Modified Julian Date UT1
 *
 * @return Greenwich Hour Angle matrix
 *
 */
//------------------------------------------------------------------------------
Matrix GHAMatrix(double Mjd_UT1) {
    return R_z( gast(Mjd_UT1) );
}
