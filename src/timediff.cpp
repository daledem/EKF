//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#include "../include/timediff.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// void timediff(double& UT1_TAI,double& UTC_GPS,double& UT1_GPS,double& TT_UTC,
//                  double& GPS_UTC,double UT1_UTC,double TAI_UTC)
//------------------------------------------------------------------------------
/**
 *   Time differences [s]
 *
 * @param[out] <UT1_TAI> UT1-TAI time difference [s]
 * @param[out] <UTC_GPS> UTC_GPS time difference [s]
 * @param[out] <UT1_GPS> UT1-GPS time difference [s]
 * @param[out] <TT_UTC> TT-UTC time difference [s]
 * @param[out] <GPS_UTC> GPS-UTC time difference [s]
 * @param <UT1_UTC> UT1-UTC time difference [s]
 * @param <TAI_UTC> TAI-UTC time difference [s]
 *
 */
//------------------------------------------------------------------------------
void timediff(double& UT1_TAI,double& UTC_GPS,double& UT1_GPS,double& TT_UTC,double& GPS_UTC,double UT1_UTC,double TAI_UTC) {
    double TT_TAI  = +32.184;          // TT-TAI time difference [s]

    double GPS_TAI = -19.0;            // GPS-TAI time difference [s]

    double TT_GPS  =  TT_TAI-GPS_TAI;  // TT-GPS time difference [s]

    double TAI_GPS = -GPS_TAI;         // TAI-GPS time difference [s]

    UT1_TAI = UT1_UTC-TAI_UTC;  // UT1-TAI time difference [s]

    double UTC_TAI = -TAI_UTC;         // UTC-TAI time difference [s]
  
    UTC_GPS = UTC_TAI-GPS_TAI;  // UTC_GPS time difference [s]

    UT1_GPS = UT1_TAI-GPS_TAI;  // UT1-GPS time difference [s]

    TT_UTC  = TT_TAI-UTC_TAI;   //  TT-UTC time difference [s]

    GPS_UTC = GPS_TAI-UTC_TAI;  // GPS-UTC time difference [s]

}
