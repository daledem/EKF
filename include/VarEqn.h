//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/11
//
//------------------------------------------------------------------------------
#ifndef VAREQN_H
#define VAREQN_H

#include "./Matrix.h"
#include "./SAT_Const.h"
#include "./Globals.h"
#include "./IERS.h"
#include "./timediff.h"
#include "./PrecMatrix.h"
#include "./NutMatrix.h"
#include "./PoleMatrix.h"
#include "./GHAMatrix.h"
#include "./AccelHarmonic.h"
#include "./G_AccelHarmonic.h"

Matrix VarEqn(double x,const Matrix& yPhi);

#endif //VAREQN_H
