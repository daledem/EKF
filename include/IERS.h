#ifndef IERS_H
#define IERS_H

#include <cmath>
#include "./Matrix.h"
#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// IERS(double& x_pole,double& y_pole,double& UT1_UTC,double& LOD,double& dpsi,double& deps,double& dx_pole,double& dy_pole,double& TAI_UTC,const Matrix& eop,double Mjd_UTC,char interp = 'n')
//------------------------------------------------------------------------------
/**
*   Management of IERS time and polar motion data
*
* @param[out] <x_pole> Pole coordinate [rad]
* @param[out] <y_pole> Pole coordinate [rad]
* @param[out] <UT1_UTC> UT1-UTC time difference [s]
* @param[out] <LOD> Length of day [s]
* @param[out] <dpsi>
* @param[out] <deps>
* @param[out] <dx_pole> Pole coordinate [rad]
* @param[out] <dy_pole> Pole coordinate [rad]
* @param[out] <TAI_UTC> TAI-UTC time difference [s]
* @param <eop> Earth orientation parameters
* @param <Mjd_UTC> Modified Julian Date UTC
* @param <interp>
*
*/
//------------------------------------------------------------------------------

void IERS(double& x_pole,double& y_pole,double& UT1_UTC,double& LOD,double& dpsi,double& deps,double& dx_pole,double& dy_pole,double& TAI_UTC,const Matrix& eop,double Mjd_UTC,char interp = 'n');

#endif //IERS_H
