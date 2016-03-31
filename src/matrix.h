#ifndef MATRIX_H
#define MATRIX_H

#include <math.h>
#include <string.h>
#include "vector.h"
#include <iostream>
#include <iomanip>

class matrix22 {
  private:
  protected:
  public:
    double a, b, c, d;
    matrix22();
    matrix22(double a, double b, double c, double d);
    matrix22 operator*(const matrix22& m);
    vector2 operator*(const vector2& v);
    matrix22 operator+(const matrix22& m);
    matrix22 operator-(const matrix22& m);
    matrix22 operator*(const double s);
    matrix22& operator=(const matrix22& m);
    double det();
    matrix22 inv();
    double trace();
    void eigenvalues(double& lambda0, double& lambda1);
};

double* gaussJordanElimination(const double m_[], int n, bool& singular);

class cMatrix {
  private:
  protected:
  public:
	double *A;
	int m, n;

	cMatrix(const int m, const int n);
	cMatrix(const int m, const int n, const double A[]);

	~cMatrix();

	cMatrix& set(const double A[]);
	cMatrix& eye();
	cMatrix& zero();
	cMatrix transpose();
	cMatrix  operator*(const cMatrix& B);
	cMatrix  operator*(const double s);
	cMatrix  operator+(const cMatrix& B);
	cMatrix  operator-(const cMatrix& B);
	cMatrix& operator=(const cMatrix& B);
	vector3 mult(const vector3& v);
	vector3 mult_(const vector3& v);
	void householderDecomposition(cMatrix& Q, cMatrix& R);
	void householderBidiagonalization(cMatrix& Q, cMatrix& R, cMatrix& S);
	void output();
};

#endif