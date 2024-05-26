//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#ifndef PROYECTO_ELEMENTS_H
#define PROYECTO_ELEMENTS_H

#include <cmath>
#include "./Matrix.h"
#include "./SAT_Const.h"
#include "./mod.h"

void elements(double& p, double& a, double& e, double& i, double& Omega, double& omega, double& M, const Matrix& y);

#endif //PROYECTO_ELEMENTS_H
