//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/13
//
//------------------------------------------------------------------------------
#ifndef MEASUPDATE_H
#define MEASUPDATE_H

#include "./Matrix.h"

void MeasUpdate(Matrix& K,Matrix& x, double z,double g,double s,const Matrix& G,Matrix& P,int n);

#endif //MEASUPDATE_H
