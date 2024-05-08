#include <assert.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>

#include "./include/Matrix.h"
#include "./include/R_x.h"
#include "./include/R_y.h"
#include "./include/R_z.h"
#include "./include/sign_.h"
#include "./include/Globals.h"
#include "./include/unit.h"
#include "./include/timediff.h"
#include "./include/Position.h"
#include "./include/Mjday_TDB.h"
#include "./include/Mjday.h"
#include "./include/MeanObliquity.h"
#include "./include/Legendre.h"
#include "./include/Geodetic.h"
#include "./include/Frac.h"
#include "./include/EccAnom.h"
#include "./include/NutAngles.h"
#include "./include/AccelPointMass.h"
#include "./include/AzElPa.h"
#include "./include/Cheb3D.h"
#include "./include/IERS.h"
#include "./include/gmst.h"

#define TOL_ 10e-14

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;


int proMat_01()
{
    double v1[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {1.0, 0.0, 0.0, 1.0};
    Matrix m1(2, 2, v1, 4);
    Matrix m2(2, 2, v2, 4);
    Matrix sol(2, 2);
    
    sol = m1 * m2;

    m1.print();
    m2.print();
    sol.print();

    _assert(sol(1,1) == m1(1,1) && sol(1,2) == m1(1,2) && sol(2,1) == m1(2,1) && sol(2,2) == m1(2,2));
    
    return 0;
}

int R_x_01(){
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_x(alpha);

    _assert(fabs(sol(1,1) - 1) < TOL_ && fabs(sol(1,2)) < TOL_ && fabs(sol(1,3)) < TOL_);
    _assert(fabs(sol(2,1)) < TOL_ && fabs(sol(2,2) - 0.540302305868140 ) < TOL_ && fabs(sol(2,3) - 0.841470984807897 ) < TOL_);
    _assert(fabs(sol(3,1)) < TOL_ && fabs(sol(3,2) + 0.841470984807897 ) < TOL_ && fabs(sol(3,3) - 0.540302305868140 ) < TOL_);

    return 0;
}

int R_y_01(){
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_y(alpha);

    _assert(fabs(sol(1,1) - 0.540302305868140 ) < TOL_ && fabs(sol(1,2)) < TOL_ && fabs(sol(1,3) + 0.841470984807897 ) < TOL_);
    _assert(fabs(sol(2,1)) < TOL_ && fabs(sol(2,2) - 1.0 ) < TOL_ && fabs(sol(2,3)) < TOL_);
    _assert(fabs(sol(3,1) - 0.841470984807897 ) < TOL_ && fabs(sol(3,2)) < TOL_ && fabs(sol(3,3) - 0.540302305868140 ) < TOL_);

    return 0;
}

int R_z_01(){
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_z(alpha);

    _assert(fabs(sol(1,1) - 0.540302305868140 ) < TOL_ && fabs(sol(1,2) - 0.841470984807897 ) < TOL_ && fabs(sol(1,3)) < TOL_);
    _assert(fabs(sol(2,1) + 0.841470984807897 ) < TOL_ && fabs(sol(2,2) - 0.540302305868140 ) < TOL_ && fabs(sol(2,3)) < TOL_);
    _assert(fabs(sol(3,1)) < TOL_ && fabs(sol(3,2)) < TOL_ && fabs(sol(3,3) - 1 ) < TOL_);

    return 0;
}

int sign__01(){
    _assert(fabs(sign_(5,-1) + 5) < TOL_);
    _assert(fabs(sign_(5,1) - 5) < TOL_);
    _assert(fabs(sign_(123,-1) + 123) < TOL_);
    _assert(fabs(sign_(123,1) - 123) < TOL_);

    return 0;
}

int unit_01(){
    double vec [] = {1,2,3};
    Matrix entrada = Matrix(1,3,vec,3);

    Matrix sol(1,3);

    sol = unit(entrada);

    double vecRes [] = {0.267261241912424,0.534522483824849,0.801783725737273};
    Matrix res = Matrix(1,3,vecRes,3);

    _assert(sol.equals(res,TOL_));

    return 0;
}

int unit_02(){
    double vec[] = {1,2,3,5};
    Matrix entrada = Matrix(1,4,vec,4);

    Matrix sol(1,3);
    sol = unit(entrada);

    double vecRes [] = {0.160128153805087,0.320256307610174,0.480384461415261};
    Matrix res = Matrix(1,3,vecRes,3);

    _assert(sol.equals(res,TOL_));

    return 0;
}

int timediff_01(){
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

    timediff(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC,5,6);

    _assert(fabs(UT1_TAI + 1) < TOL_);
    _assert(fabs(UTC_GPS - 13) < TOL_);
    _assert(fabs(UT1_GPS - 18) < TOL_);
    _assert(fabs(TT_UTC - 38.183999999999997) < TOL_);
    _assert(fabs(GPS_UTC + 13) < TOL_);

    return 0;
}

int position_01(){
    Matrix sol(1,3);
    sol = Position(1,2,3);
    double matriz[] = {-1.438078785611559*1.0e06,-2.239675009373783*1.0e06,5.776810445003163*1.0e06};
    Matrix expected(1,3,matriz,3);

    _assert(sol.equals(expected,TOL_));


    return 0;
}

int Mjday_TDB_01(){
    double sol;
    sol = Mjday_TDB(123);

    _assert(fabs(sol - 1.230000000183996e+02) < TOL_);

    return 0;
}

int Mjday_01(){
    double sol;
    sol = Mjday(1,2,3);

    _assert(fabs(sol + 678557) < TOL_);

    return 0;
}

int Mjday_02(){
    double sol;
    sol = Mjday(1,2,3,4,5,6);

    _assert(fabs(sol + 6.785568297916667e+05) < TOL_);

    return 0;
}

int MeanObliquity_01(){
    double sol;
    sol = MeanObliquity(123);

    _assert(fabs(sol - 0.409412306065671) < TOL_);

    return 0;
}

int Legendre_01(){
    Matrix pnm(2,2);
    Matrix dpnm(2,2);
    Legendre(pnm,dpnm,1,1,5);

    double matriz0 [] = {1.0,0.0,-1.66090556432769,0.49131731740833};
    Matrix expected0(2,2,matriz0,4);

    double matriz1 [] = {0.0,0.0,0.49131731740833,1.66090556432769};
    Matrix expected1(2,2,matriz1,4);

    _assert(pnm.equals(expected0,TOL_));
    _assert(dpnm.equals(expected1,TOL_));

    return 0;
}


int Geodetic_01(){
    double lon,lat,h;
    double vr[]= {1,2,3};
    Matrix r(1,3,vr,3);
    Geodetic(lon,lat,h,r);

    _assert(fabs(lon - 1.107148717794090) < TOL_);
    _assert(fabs(lat - 1.570744136243924) < TOL_);
    _assert(fabs(h + 6.356748616533795e+06) < TOL_);

    return 0;
}

int Frac_01(){
    double sol;
    sol = Frac(1.26);

    _assert(fabs(sol - 0.26) < TOL_);

    return 0;
}

int EccAnom_01(){
    double sol;
    sol = EccAnom(12,0.5);

    _assert(fabs(sol - 5.300942369055246) < TOL_);

    return 0;
}


int NutAngles_01(){
    double dpsi,deps;
    NutAngles(dpsi,deps,12);

    _assert(fabs(dpsi - 3.145583548136575e-05) < TOL_);
    _assert(fabs(deps - 3.844741370508873e-05) < TOL_);

    return 0;
}


int AccelPointMass_01(){
    double v1[] = {1,2,3};
    Matrix r(1,3,v1,3);
    double v2[] = {4,5,6};
    Matrix s(1,3,v2,3);

    Matrix sol(1,3);
    sol = AccelPointMass(r,s,7);

    double matriz [] = {0.108243193501550,0.097883141096128,0.087523088690707};
    Matrix expected(1,3,matriz,3);

    _assert(sol.equals(expected,TOL_));


    return 0;
}


int AzElPa_01(){
    double Az,El;
    Matrix dAds(1,3);
    Matrix dEds(1,3);
    double v[] = {1,2,3};
    Matrix s(1,3,v,3);

    AzElPa(Az,El,dAds,dEds,s);

    double v1[] = {0.400000000000000,-0.200000000000000,0};
    Matrix expected1(1,3,v1,3);
    double v2[] = {-0.095831484749991,-0.191662969499982,0.159719141249985};
    Matrix expected2(1,3,v2,3);

    _assert(fabs(Az - 0.463647609000806) < TOL_);
    _assert(fabs(El - 0.930274014115472) < TOL_);
    _assert(dAds.equals(expected1,TOL_));
    _assert(dEds.equals(expected2,TOL_));

    return 0;
}


int Cheb3D_01() {
    double vx[] = {1,2,3};
    Matrix Cx(1,3,vx,3);

    double vy[] = {4,5,6};
    Matrix Cy(1,3,vy,3);

    double vz[] = {7,8,9};
    Matrix Cz(1,3,vz,3);

    Matrix sol(1,3);

    sol = Cheb3D(3,3,2,4,Cx,Cy,Cz);

    double v[] = {-2,-2,-2};
    Matrix expected(1,3,v,3);

    _assert(sol.equals(expected,TOL_));

    return 0;
}

int IERS_01() {
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,*Global::eopdata,37666.00000000000000,'l');

    _assert(fabs(x_pole + 7.708537529641623e-08) < TOL_);
    _assert(fabs(y_pole - 1.037986091255516e-06) < TOL_);
    _assert(fabs(UT1_UTC - 0.032054700000000) < TOL_);
    _assert(fabs(LOD - 0.001669000000000) < TOL_);
    _assert(fabs(dpsi - 3.101789450370700e-07) < TOL_);
    _assert(fabs(deps - 3.049478054178981e-08) < TOL_);
    _assert(fabs(dx_pole) < TOL_);
    _assert(fabs(dy_pole) < TOL_);
    _assert(fabs(TAI_UTC - 2) < TOL_);

    return 0;
}


int IERS_02() {
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,*Global::eopdata,37666);

    _assert(fabs(x_pole + 7.708537529641623e-08) < TOL_);
    _assert(fabs(y_pole - 1.037986091255516e-06) < TOL_);
    _assert(fabs(UT1_UTC - 0.032054700000000) < TOL_);
    _assert(fabs(LOD - 0.001669000000000) < TOL_);
    _assert(fabs(dpsi - 3.101789450370700e-07) < TOL_);
    _assert(fabs(deps - 3.049478054178981e-08) < TOL_);
    _assert(fabs(dx_pole) < TOL_);
    _assert(fabs(dy_pole) < TOL_);
    _assert(fabs(TAI_UTC - 2) < TOL_);

    return 0;
}


int gmst_01() {
    double sol;
    sol = gmst(37666);

    _assert(fabs(sol - 1.765465139678228) < TOL_);

    return 0;
}


int all_tests()
{
    _verify(proMat_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(sign__01);
    _verify(unit_01);
    _verify(unit_02);
    _verify(timediff_01);
    _verify(position_01);
    _verify(Mjday_TDB_01);
    _verify(Mjday_01);
    _verify(Mjday_02);
    _verify(MeanObliquity_01);
    _verify(Legendre_01);
    _verify(Geodetic_01);
    _verify(Frac_01);
    _verify(EccAnom_01);
    _verify(NutAngles_01);
    _verify(AccelPointMass_01);
    _verify(AzElPa_01);
    _verify(Cheb3D_01);
    _verify(IERS_01);
    _verify(IERS_02);
    _verify(gmst_01);
 
    return 0;
}


int main()
{
    Global::eop19620101(6);
    Global::eopdata->print();

    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}

