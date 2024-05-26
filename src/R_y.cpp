//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/18
//
//------------------------------------------------------------------------------
#include "../include/R_y.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix R_y(double angle)
//------------------------------------------------------------------------------
/**
 * Calculates rotation matrix for the y direction.
 *
 * @param <angle> angle of rotation [rad]
 *
 * @return vector result
 */
//------------------------------------------------------------------------------
Matrix R_y(double angle)
{

    double C = cos(angle);
    double S = sin(angle);
    Matrix rotmat(3,3);

    rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
    rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
    rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;

    return rotmat;
}