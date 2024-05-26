//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/09
//
//------------------------------------------------------------------------------
#include "../include/PoleMatrix.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix NutMatrix(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 *   Transformation from pseudo Earth-fixed to Earth-fixed coordinates
 *      for a given date
 *
 * @param <xp> Pole coordinte
 * @param <yp> Pole coordinte
 *
 * @return Pole matrix
 *
 */
//------------------------------------------------------------------------------
Matrix PoleMatrix(double xp, double yp) {
    return R_y(-xp) * R_x(-yp);
}
