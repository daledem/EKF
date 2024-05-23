#ifndef HGIBBS_H
#define HGIBBS_H

#include <cmath>
#include <string>
#include "./Matrix.h"
#include "./SAT_Const.h"
#include "./unit.h"
#include "./angl.h"

void hgibbs(Matrix& v2,double& theta,double& theta1,double& copa,std::string& error,const Matrix& r1,const Matrix& r2,const Matrix& r3,double Mjd1,double Mjd2,double Mjd3);

#endif //HGIBBS_H
