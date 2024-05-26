//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/11
//
//------------------------------------------------------------------------------
#ifndef TIMEUPDATE_H
#define TIMEUPDATE_H

#include "./Matrix.h"

Matrix TimeUpdate(const Matrix& P,const Matrix& Phi, double Qdt = 0.0);

#endif //TIMEUPDATE_H
