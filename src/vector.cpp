#include "../include/vector.h"

// Función que calcula la norma de un vector
// In: v vector (double)
// In: n dimensión (int)
// Out: devuelve la norma de v
double norm(double v[], int n)
{
    double suma = 0.0;
    int i;

    if(n <= 0)
        throw "Empty vector";

    for(i = 0; i < n; ++i)
        suma += v[i]*v[i];

    return(sqrt(suma));
}


void operation(double result[], int &tam, double operando1[], int tam1, double operando2, double (*func)(double &i, double &j)) {
    if((tam1 <= 0))
        throw "Empty vector";

    for(int i = 0;i < tam1;i++)
        result[i] =  func(operando1[i],operando2);

    tam = tam1;
}


void operation(double result[], int &tam, double operando1[], int tam1, double operando2[], int tam2, double (*func)(double &i, double &j)) {
    if((tam1 <= 0) || (tam2 <= 0) || (tam1 != tam2))
        throw "Empty vector or different dimensions";

    for(int i = 0;i < tam1;i++)
        result[i] =  func(operando1[i],operando2[i]);

    tam = tam1;
}


double sum(double &i, double &j) {
    return i + j;
}


double mul(double &i, double &j) {
    return i*j;
}


double div(double &i, double &j) {
    return i/j;
}
