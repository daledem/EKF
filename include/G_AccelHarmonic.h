//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/09
//
//------------------------------------------------------------------------------
#ifndef G_ACCELHARMONIC_H
#define G_ACCELHARMONIC_H

#include "./Matrix.h"
#include "./AccelHarmonic.h"

Matrix G_AccelHarmonic(const Matrix& r,const Matrix& U,int n_max,int m_max);

#endif //G_ACCELHARMONIC_H
