//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#ifndef NUTANGLES_H
#define NUTANGLES_H

#include "./Matrix.h"
#include "./SAT_Const.h"
#include <cmath>

void NutAngles(double& dpsi,double& deps,double Mjd_TT);

#endif //NUTANGLES_H
