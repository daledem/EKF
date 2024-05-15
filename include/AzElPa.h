#ifndef AZELPA_H
#define AZELPA_H

#include <cmath>
#include "./Matrix.h"
#include "./SAT_Const.h"

void AzElPa(double& Az,double& El,Matrix& dAds,Matrix& dEds,const Matrix& s);

#endif //AZELPA_H
