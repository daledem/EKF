#ifndef ACCEL_H
#define ACCEL_H

#include "./Matrix.h"
#include "./SAT_Const.h"
#include "./Globals.h"
#include "./IERS.h"
#include "./timediff.h"
#include "./PrecMatrix.h"
#include "./NutMatrix.h"
#include "./PoleMatrix.h"
#include "./GHAMatrix.h"
#include "./Mjday_TDB.h"
#include "./JPL_Eph_DE430.h"
#include "./AccelHarmonic.h"
#include "./AccelPointMass.h"

Matrix Accel(double x,const Matrix& Y);

#endif //ACCEL_H
