#ifndef _MATRIX_
#define _MATRIX_

class Matrix
{
    public:
        Matrix(int fil, int col);
        Matrix(int fil, int col, double v[], int n);
        Matrix(const Matrix& m);
        ~Matrix();
 
        Matrix& operator=(const Matrix& matrix2);
        Matrix  operator+(const Matrix& matrix2);
        Matrix  operator+(double sumando);
        friend Matrix operator+(const double& escalar,const Matrix& matrix);
        Matrix  operator-(const Matrix& matrix2) const;
        Matrix  operator-();
        Matrix  operator*(const Matrix& matrix2) const;
        Matrix  operator*(double multiplicador) const;
        friend Matrix operator*(const double& escalar,const Matrix& matrix);
        Matrix  operator/(double divisor);
        double& operator()(const int i, const int j) const;

        Matrix trans() const;
        Matrix append(Matrix& matrix2);
        Matrix join(Matrix& matrix2);
        int getFil() const;
        int getCol() const;
        Matrix getFilaByIndex(int fil,int inicio = 1) const;
        Matrix getFilaByIndex(int fil,int inicio,int fin) const;
        Matrix getColumnaByIndex(int col,int inicio = 1) const;
        Matrix getColumnaByIndex(int col,int inicio,int fin) const;
        bool equals(const Matrix& matrix2,const double TOL) const;
        void print();

        static double dot(const Matrix& matrix1,const Matrix& matrix2);
        static double norm(const Matrix& matrix);
 
    private:
        void initMatrix();
 
    private:
        int fil;
        int col;
        double **matrix;
};

#endif
