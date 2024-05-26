//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#ifndef POSITION_H
#define POSITION_H

#include "./Matrix.h"
#include "./SAT_Const.h"
#include <cmath>

Matrix Position(double lon, double lat, double h);

#endif //POSITION_H
