#ifndef TIMEDIFF_H
#define TIMEDIFF_H

//------------------------------------------------------------------------------
// timediff(double& UT1_TAI,double& UTC_GPS,double& UT1_GPS,double& TT_UTC,double& GPS_UTC,double UT1_UTC,double TAI_UTC)
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

void timediff(double& UT1_TAI,double& UTC_GPS,double& UT1_GPS,double& TT_UTC,double& GPS_UTC,double UT1_UTC,double TAI_UTC);

#endif //TIMEDIFF_H
