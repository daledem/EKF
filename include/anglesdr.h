#ifndef ANGLESDR_H
#define ANGLESDR_H

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
#include "./gibbs.h"
#include "./doubler.h"

void anglesdr(Matrix& r2,Matrix& v2,double az1,double az2,double az3,double el1,double el2,double el3,double Mjd1,double Mjd2,double Mjd3,Matrix rsite1,Matrix rsite2,Matrix rsite3);

#endif //ANGLESDR_H
