//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/09
//
//------------------------------------------------------------------------------
#ifndef PRECMATRIX_H
#define PRECMATRIX_H

#include "./Matrix.h"
#include "./R_z.h"
#include "./R_y.h"
#include "./SAT_Const.h"

Matrix PrecMatrix(double Mjd_1,double Mjd_2);

#endif //PRECMATRIX_H
