#ifndef GEODETIC_H
#define GEODETIC_H

#include <cmath>
#include "./Matrix.h"
#include "./SAT_Const.h"

void Geodetic(double& lon, double& lat, double& h, const Matrix& r);

#endif //GEODETIC_H
