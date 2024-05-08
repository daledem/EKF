#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#define TOL_ 10e-14

Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}
 
Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();
 
    int k = 0;
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}
 
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}
 
Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}
 
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];
 
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}
 
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];
 
    return *this;
}
 
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, col);
 
    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
 
    return result;
}


Matrix Matrix::operator*(double multiplicador)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] * multiplicador;

    return result;
}


Matrix operator*(const double &escalar, const Matrix &matrix) {
    Matrix result(matrix.fil, matrix.col);

    for (int i = 0; i < matrix.fil; i++)
        for (int j = 0; j < matrix.col; j++)
            result.matrix[i][j] = matrix.matrix[i][j] * escalar;

    return result;
}


Matrix Matrix::operator/(double divisor)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] / divisor;

    return result;
}

 
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}


Matrix Matrix::getFil(int fil) const {
    Matrix result(1,col);

    for (int i = 0; i < col; i++)
        result.matrix[0][i] = matrix[fil - 1][i];

    return result;
}


Matrix Matrix::getCol(int col) const {
    Matrix result(fil, 1);

    for (int i = 0; i < fil; i++)
        result.matrix[i][0] = matrix[i][col - 1];

    return result;
}


bool Matrix::equals(const Matrix &matrix2,const double TOL) const {

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            if(fabs(matrix[i][j] - matrix2.matrix[i][j]) > TOL)
                return false;
    return true;
}


void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


int Matrix::find(const double objective, const Matrix &matrix,const int fil) {
    int i = 0;
    int j = 1;

    while(j < matrix.col && i == 0){
        if(fabs(objective-matrix(fil,j)) < 10e-14) {
            i = j;
        }
        j++;
    }

    return i;
}


double Matrix::dot(const Matrix &matrix1, const Matrix &matrix2) {
    if(matrix1.col != matrix2.col)
        throw "Different dimensions";

    double sum = 0.0;

    for (int j = 1; j <= matrix1.col; j++){
        sum += matrix1(1,j)*matrix2(1,j);
    }

    return sum;
}


double Matrix::norm(const Matrix& matrix) {
    double sum = 0.0;

    for (int j = 1; j <= matrix.col; j++){
        sum += matrix(1,j)*matrix(1,j);
    }

    return sqrt(sum);
}
