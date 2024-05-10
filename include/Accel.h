#ifndef ACCEL_H
#define ACCEL_H

#include "./Matrix.h"
#include "./SAT_Const.h"
#include "./Globals.h"
#include "./IERS.h"
#include "./timediff.h"
#include "./PrecMatrix.h"


Matrix Accel(double x,const Matrix& Y);

#endif //ACCEL_H
