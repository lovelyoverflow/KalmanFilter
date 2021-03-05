#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>

class Matrix 
{
public:
	Matrix();
	Matrix(std::vector<std::vector<double>> m);
	bool Zero(int dimX, int dimY);
	bool Identity(int dim);
	void Show();
	Matrix operator+(Matrix pMat);
	Matrix operator-(Matrix pMat);
	Matrix operator*(Matrix pMat);
	Matrix Transpose();
	Matrix Cholesky(double ztol);
	Matrix CholeskyInverse();
	Matrix Inverse();
	double GetMean();
	std::vector<std::vector<double>> mat;
private:
	int dimX;
	int dimY;
};

#endif