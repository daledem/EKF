//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#ifndef ANGL_H
#define ANGL_H

#include <cmath>
#include "./Matrix.h"
#include "./sign.h"

double angl(const Matrix& vec1, const Matrix& vec2);

#endif //ANGL_H
