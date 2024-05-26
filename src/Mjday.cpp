//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/28
//
//------------------------------------------------------------------------------
#include "../include/Mjday.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// double Mjday(double yr, double mon, double day, double hr = 0, double min  = 0, double sec = 0)
//------------------------------------------------------------------------------
/**
 *   Computes the modified julian date of a certain date
 *
 * @param <yr> year
 * @param <mon> month
 * @param <day> day
 * @param <hr> universal time hour
 * @param <min> universal time min
 * @param <sec> universal time sec
 *
 * @return Modified julian date
 *
 */
//------------------------------------------------------------------------------
double Mjday(double yr, double mon, double day, double hr, double min, double sec) {

    double jd = 367.0 * yr
        - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )
        + floor( 275 * mon / 9.0 )
        + day + 1721013.5
        + ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;

    return jd-2400000.5;
}

