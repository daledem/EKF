//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <cmath>
#include "./Matrix.h"

void Legendre(Matrix& pnm,Matrix& dpnm,int n, int m, double fi);

#endif //LEGENDRE_H
