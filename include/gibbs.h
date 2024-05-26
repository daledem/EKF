//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#ifndef GIBBS_H
#define GIBBS_H

#include <cmath>
#include <string>
#include "./Matrix.h"
#include "./SAT_Const.h"
#include "./unit.h"
#include "./angl.h"

void gibbs(Matrix& v2,double& theta,double& theta1,double& copa,std::string& error,const Matrix& r1,const Matrix& r2,const Matrix& r3);

#endif //GIBBS_H
