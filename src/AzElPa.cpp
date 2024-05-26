//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/02
//
//------------------------------------------------------------------------------
#include "../include/AzElPa.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// void AzElPa(double& Az,double& El,Matrix& dAds,Matrix& dEds,const Matrix& s)
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
void AzElPa(double &Az, double &El, Matrix &dAds, Matrix &dEds, const Matrix &s) {
    double pi2 = 2.0*Const::pi;

    double rho = sqrt(s(1,1)*s(1,1)+s(2,1)*s(2,1));

    // Angles
    Az = atan2(s(1,1),s(2,1));

    if (Az<0.0)
        Az = Az+pi2;

    El = atan ( s(3,1) / rho );

    // Partials
    double vdAds[] = { s(2,1)/(rho*rho), -s(1,1)/(rho*rho), 0.0 };
    double vdEds[] = { -s(1,1)*s(3,1)/rho, -s(2,1)*s(3,1)/rho , rho };

    dAds = Matrix(1,3,vdAds,3);
    dEds = Matrix(1,3,vdEds,3) / Matrix::dot(s,s);
}
