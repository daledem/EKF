//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#ifndef ANGLESG_H
#define ANGLESG_H

#include <cmath>
#include "./Matrix.h"
#include "./SAT_Const.h"
#include "./Geodetic.h"
#include "./LTC.h"
#include "./IERS.h"
#include "./timediff.h"
#include "./PrecMatrix.h"
#include "./NutMatrix.h"
#include "./PoleMatrix.h"
#include "./GHAMatrix.h"
#include "./rpoly.h"
#include "./gibbs.h"
#include "./hgibbs.h"
#include "./elements.h"

void anglesg(Matrix& r2,Matrix& v2,double az1,double az2,double az3,double el1,double el2,double el3,double Mjd1,double Mjd2,double Mjd3,Matrix Rs1,Matrix Rs2,Matrix Rs3);

#endif //ANGLESG_H
