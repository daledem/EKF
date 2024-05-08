#include "../include/MeanObliquity.h"

double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT-Const::MJD_J2000)/36525;

    return Const::Rad *( 84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600 );
}

