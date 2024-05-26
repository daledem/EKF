//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#include "../include/Geodetic.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// void Geodetic(double& lon, double& lat, double& h, const Matrix& r)
//------------------------------------------------------------------------------
/**
 *   geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 *   from given position vector (r [m])
 *
 * @param <lon> Longitude [rad]
 * @param <lat> latitude [rad]
 * @param <h> altitude [m]
 * @param <r> position vector [m]
 *
 * @return Geodetic coordinates
 *
 */
//------------------------------------------------------------------------------
void Geodetic(double& lon, double& lat, double& h,const Matrix& r) {
    double R_equ = Const::R_Earth;
    double f     = Const::f_Earth;

    double epsRequ = Const::eps*R_equ;         // Convergence criterion
    double e2      = f*(2.0-f);         // Square of eccentricity

    double X = r(1,1);                    // Cartesian coordinates
    double Y = r(1,2);
    double Z = r(1,3);
    double rho2 = X*X + Y*Y;            // Square of distance from z-axis

     // Check validity of input data
    if (Matrix::norm(r)==0.0) {

        lon = 0.0;
        lat = 0.0;
        h   = -Const::R_Earth;
        throw "invalid input in Geodetic constructor\n";
    }


     // Iteration
    double dZ = e2*Z;
    double ZdZ,Nh,SinPhi,N,dZ_new;

    while(1) {
        ZdZ    =  Z + dZ;
        Nh     =  sqrt( rho2 + ZdZ*ZdZ );
        SinPhi =  ZdZ / Nh;                     // Sine of geodetic latitude
        N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        dZ_new =  N*e2*SinPhi;
        if ( fabs(dZ-dZ_new) < epsRequ )
            break;

        dZ = dZ_new;
    }

     // Longitude, latitude, altitude
    lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;
}
