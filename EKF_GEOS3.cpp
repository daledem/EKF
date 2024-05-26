#include <iostream>
#include <math.h>
#include <string.h>

#include "./include/TimeUpdate.h"
#include "./include/VarEqn.h"
#include "./include/Position.h"
#include "./include/Mjday.h"
#include "./include/AzElPa.h"
#include "./include/Accel.h"
#include "./include/LTC.h"
#include "./include/DEInteg.h"
#include "./include/MeasUpdate.h"
#include "./include/anglesg.h"
#include "./include/anglesdr.h"

using namespace std;

Matrix Snm(181,181);
Matrix Cnm(181,181);
Matrix PC(2285,1020);
Matrix eopdata(13,21413);
struct auxParam AuxParam;

int main() {
     double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,
         t_old,n_eqn,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,Mjd_TT,
         Mjd_UT1,theta,Azim,Elev,Dist;
     int nobs;
     Matrix K;

     // Cargamos PC con los datos de DE430Coeff
     FILE *fid = fopen("data/M_tab.txt","r");

     if(fid == NULL){
         printf("Fail to open DE430Coeff.txt");
         exit(EXIT_FAILURE);
     }

     for(int i = 1; i <= 2285; i++){
         for(int j = 1; j <= 1020; j++){
             fscanf(fid,"%lf",&(PC(i,j)));

         }
     }

     fclose(fid);


     // Cargamos Snm y Cnm con los datos de GGM03S
     fid = fopen("data/GGM03S.txt","r");
     int aux1,aux2;
     double aux3,aux4;
     if(fid == NULL){
         printf("Fail to open GGM03S.txt");
         exit(EXIT_FAILURE);
     }

     for(int i = 0; i <= 180; i++){
         for(int j = 0; j <= i; j++){
             fscanf(fid,"%d %d %lf %lf %lf %lf",&aux1,
                    &aux2,&(Cnm(i+1,j+1)),&(Snm(i+1,j+1)),&aux3,
                    &aux4);
         }
     }

     fclose(fid);

     // Model parameters
     AuxParam.Mjd_UTC = 0;
     AuxParam.n = 0;
     AuxParam.m = 0;

     // Cargamos eopdata con los datos de eop19620101
     fid = fopen("data/eop19620101.txt","r");

     if(fid == NULL){
         printf("Fail to open eop19620101.txt");
         exit(EXIT_FAILURE);
     }

     for(int i = 1; i<=21413; i++){
         fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(eopdata(1,i)),
                &(eopdata(2,i)),&(eopdata(3,i)),&(eopdata(4,i)),&(eopdata(5,i)),
                &(eopdata(6,i)),&(eopdata(7,i)),&(eopdata(8,i)),&(eopdata(9,i)),
                &(eopdata(10,i)),&(eopdata(11,i)),&(eopdata(12,i)),&(eopdata(13,i)));
     }

     fclose(fid);

     nobs = 46;
     Matrix obs(nobs,4);

     fid = fopen("data/GEOS3.txt","r");

     if(fid == NULL) {
         printf("Fail to open GEOS3.txt");
         exit(EXIT_FAILURE);
     }

     char tline[56];
     char yy[56],mm[56],dd[56],hh[56],min[56],ss[56],azz[56],ell[56],dist[56];
     int i = 1;
     while (fgets(tline, sizeof(tline), fid) && i <= nobs) {

         strncpy(yy,tline,4);
         int Y = atoi(yy);

         strncpy(mm,tline + 5,7);
         int M = atoi(mm);

         strncpy(dd,tline + 8,10);
         int D = atoi(dd);

         strncpy(hh,tline + 12,14);
         int h = atoi(hh);

         strncpy(min,tline + 15,17);
         int m = atoi(min);

         strncpy(ss,tline + 18,24);
         double s = atof(ss);

         strncpy(azz,tline + 25,33);
         double az = atof(azz);

         strncpy(ell,tline + 35,42);
         double el = atof(ell);

         strncpy(dist,tline + 44,54);
         double Dist = atof(dist);

         obs(i,1) = Mjday(Y,M,D,h,m,s);
         obs(i,2) = Const::Rad*az;
         obs(i,3) = Const::Rad*el;
         obs(i,4) = 1e3*Dist;
         i++;
     }

     fclose(fid);

     double sigma_range = 92.5;          // [m]
     double sigma_az = 0.0224*Const::Rad; // [rad]
     double sigma_el = 0.0139*Const::Rad; // [rad]

     // Kaena Point station
     double lat = Const::Rad*21.5748;     // [rad]
     double lon = Const::Rad*(-158.2706); // [rad]
     double alt = 300.20;                // [m]

     Matrix Rs = Position(lon, lat, alt).trans();

     double Mjd1 = obs(1,1);
     double Mjd2 = obs(9,1);
     double Mjd3 = obs(18,1);

    Matrix r2(3,1);
    Matrix v2(3,1);
    anglesg(r2,v2,obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),
                 Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
    //anglesdr(r2,v2,obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),
    //             Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

     Matrix Y0_apr = r2.join(v2);

     double Mjd0 = Mjday(1995,1,29,02,38,0);

     double Mjd_UTC = obs(9,1);

     AuxParam.Mjd_UTC = Mjd_UTC;
     AuxParam.n      = 20;
     AuxParam.m      = 20;
     AuxParam.sun     = 1;
     AuxParam.moon    = 1;
     AuxParam.planets = 1;

     n_eqn  = 6;

     Matrix Y = DEInteg(Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);


     Matrix P(6,6);

     for (i = 1; i<= 3; i++)
         P(i,i)=1e8;

     for (i = 4; i<= 6; i++)
         P(i,i)=1e3;


     Matrix LT = LTC(lon,lat);

     Matrix yPhi(42,1);
     Matrix Phi(6,6);

     // Measurement loop
     double t = 0;

     for (i = 1; i <= nobs; i++) {
         // Previous step
         t_old = t;
         Matrix Y_old = Y;

         // Time increment and propagation
         Mjd_UTC = obs(i,1);                       // Modified Julian Date
         t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]

         IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,Mjd_UTC,'l');
         timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
         Mjd_TT = Mjd_UTC + TT_UTC/86400;
         Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
         AuxParam.Mjd_UTC = Mjd_UTC;
         AuxParam.Mjd_TT = Mjd_TT;

         for (int ii = 1 ; ii <= 6; ii++) {
             yPhi(ii,1) = Y_old(ii,1);
             for(int j = 1 ; j <= 6; j++) {
                 if (ii==j)
                     yPhi(6*j+ii,1) = 1;
                 else
                     yPhi(6*j+ii,1) = 0;
             }
         }

         yPhi = DEInteg (VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);

         // Extract state transition matrices
         for (int j = 1 ; j <= 6; j++)
             for (int p = 1; p <= 6; p++)
                 Phi(p,j) = yPhi(6*j+p,1);

         Y = DEInteg (Accel,0,t-t_old,1e-13,1e-6,6,Y_old);

         // Topocentric coordinates
         theta = gmst(Mjd_UT1);                    // Earth rotation
         Matrix U = R_z(theta);
         Matrix r = Y.getColumnaByIndex(1,1,3);
         Matrix s = LT*(U*r-Rs);                          // Topocentric position [m]

         // Time update
         P = TimeUpdate(P, Phi);

         // Azimuth and partials
         Matrix dAds,dEds;
         AzElPa(Azim, Elev, dAds, dEds,s);     // Azimuth, Elevation
         Matrix aux(1,3);
         Matrix dAdY = (dAds*LT*U).append(aux);

         // Measurement update
         MeasUpdate (K, Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );

         // Elevation and partials
         r = Y.getColumnaByIndex(1,1,3);
         s = LT*(U*r-Rs);                          // Topocentric position [m]
         AzElPa(Azim, Elev, dAds, dEds,s);     // Azimuth, Elevation
         Matrix dEdY = (dEds*LT*U).append(aux);

         // Measurement update
         MeasUpdate (K,Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );

         // Range and partials
         r = Y.getColumnaByIndex(1,1,3);
         s = LT*(U*r-Rs);                          // Topocentric position [m]
         Dist = Matrix::norm(s);
         Matrix dDds = (s/Dist).trans();         // Range
         Matrix dDdY = (dDds*LT*U).append(aux);

         // Measurement update
         MeasUpdate (K,Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );
     }

     IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,obs(46,1),'l');
     timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
     Mjd_TT = Mjd_UTC + TT_UTC/86400;
     AuxParam.Mjd_UTC = Mjd_UTC;
     AuxParam.Mjd_TT = Mjd_TT;
     Matrix Y0 = DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

     double Y_true[] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};

     printf("\nError of Position Estimation\n");
     printf("dX %.1f [m]\n",Y0(1,1)-Y_true[0]);
     printf("dY %.1f [m]\n",Y0(2,1)-Y_true[1]);
     printf("dZ %.1f [m]\n",Y0(3,1)-Y_true[2]);
     printf("\nError of Velocity Estimation\n");
     printf("dVx %.1f [m/s]\n",Y0(4,1)-Y_true[3]);
     printf("dVy %.1f [m/s]\n",Y0(5,1)-Y_true[4]);
     printf("dVz %.1f [m/s]\n",Y0(6,1)-Y_true[5]);

     return 0;
 }
