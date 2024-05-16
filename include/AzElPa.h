#ifndef AZELPA_H
#define AZELPA_H

#include <cmath>
#include "./Matrix.h"
#include "./SAT_Const.h"

//------------------------------------------------------------------------------
// AzElPa(double& Az,double& El,Matrix& dAds,Matrix& dEds,const Matrix& s)
//------------------------------------------------------------------------------
/**
*   Computes azimuth, elevation and partials from local tangent coordinates
*
* @param[out] <Az> Azimuth [rad]
* @param[out] <El> Elevation [rad]
* @param[out] <dAds> Partials of azimuth w.r.t. s
* @param[out] <dEds> Partials of elevation w.r.t. s
* @param <s> Topocentric local tangent coordinates (East-North-Zenith frame)
*
*/
//------------------------------------------------------------------------------

void AzElPa(double& Az,double& El,Matrix& dAds,Matrix& dEds,const Matrix& s);

#endif //AZELPA_H
