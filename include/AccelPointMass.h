//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#ifndef ACCELPOINTMASS_H
#define ACCELPOINTMASS_H

#include <cmath>
#include "./Matrix.h"

Matrix AccelPointMass(const Matrix& r, Matrix& s, double GM);

#endif //ACCELPOINTMASS_H
