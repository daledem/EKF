//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#ifndef DOUBLER_H
#define DOUBLER_H

#include <cmath>
#include "./Matrix.h"
#include "./SAT_Const.h"

void doubler(Matrix& r2,Matrix& r3,double& f1,double& f2,double& q1,double& magr1,double& magr2,double& a,double& deltae32,
    double cc1,double cc2,double magrsite1,double magrsite2,double magr1in,double magr2in, const Matrix& los1,
    const Matrix& los2,const Matrix& los3,const Matrix& rsite1,const Matrix& rsite2,const Matrix& rsite3,double t1,
    double t3,char direct);

#endif //DOUBLER_H
