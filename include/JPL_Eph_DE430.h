//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/10
//
//------------------------------------------------------------------------------
#ifndef JPL_EPH_DE430_H
#define JPL_EPH_DE430_H

#include <cmath>
#include "./Matrix.h"
#include "./Cheb3D.h"

void JPL_Eph_DE430(Matrix& r_Mercury,Matrix& r_Venus,Matrix& r_Earth, Matrix& r_Mars,Matrix& r_Jupiter,Matrix& r_Saturn,
    Matrix& r_Uranus,Matrix& r_Neptune,Matrix& r_Pluto,Matrix& r_Moon,Matrix& r_Sun, double Mjd_TDB);

#endif //JPL_EPH_DE430_H
