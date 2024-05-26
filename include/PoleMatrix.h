//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/09
//
//------------------------------------------------------------------------------
#ifndef POLE_H
#define POLE_H

#include "./Matrix.h"
#include "./R_x.h"
#include "./R_y.h"

Matrix PoleMatrix(double xp,double yp);

#endif //POLE_H
