#ifndef MJDAY_H
#define MJDAY_H

#include <cmath>

//--------------------------------------------------------------------------
//  inputs:
//    year        - year                           
//    mon         - month                          
//    day         - day                            
//    hr          - universal time hour            
//    min         - universal time min             
//    sec         - universal time sec             
//
//  output:       
//    Mjd         - Modified julian date                    
//--------------------------------------------------------------------------

double Mjday(double yr, double mon, double day, double hr = 0, double min  = 0, double sec = 0);

#endif //MJDAY_H
