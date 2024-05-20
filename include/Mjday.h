#ifndef MJDAY_H
#define MJDAY_H

#include <cmath>

//------------------------------------------------------------------------------
// Mjday(double yr, double mon, double day, double hr = 0, double min  = 0, double sec = 0)
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

double Mjday(double yr, double mon, double day, double hr = 0, double min  = 0, double sec = 0);

#endif //MJDAY_H
