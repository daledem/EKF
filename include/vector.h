#ifndef PROYECTO_NORM_H
#define PROYECTO_NORM_H

#include <cmath>

double norm(double v[], int n = 3);
void operation(double result[], int &tam, double operando1[],int tam1, double operando2[],int tam2,double (*func)(double &i,double &j));
void operation(double result[], int &tam, double operando1[],int tam1, double operando2,double (*func)(double &i,double &j));
double sum(double &i,double &j);
double mul(double &i,double &j);
double div(double &i,double &j);

#endif //PROYECTO_NORM_H
