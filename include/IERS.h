//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/05
//
//------------------------------------------------------------------------------
#ifndef IERS_H
#define IERS_H

#include <cmath>
#include "./Matrix.h"
#include "./SAT_Const.h"

void IERS(double& x_pole,double& y_pole,double& UT1_UTC,double& LOD,double& dpsi,double& deps,double& dx_pole,double& dy_pole,double& TAI_UTC,const Matrix& eop,double Mjd_UTC,char interp = 'n');

#endif //IERS_H
