#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include "Matrix.h"

Matrix::Matrix()
	: mat(NULL), dimX(0), dimY(0)
{}

Matrix::Matrix(std::vector<std::vector<double>> m)
	: mat(m), dimX(m.size()), dimY(m[0].size())
{}

bool Matrix::Zero(int dimX, int dimY)
{
	if (dimX < 1 || dimY < 1)
		return false;

	this->dimX = dimX;
	this->dimY = dimY;
	
	for (int i = 0; i < dimX; i++)
	{
		std::vector<double> row;
		for (int j = 0; j < dimY; j++)
			row.push_back(0);

		this->mat.push_back(row);
	}

	return true;
}

bool Matrix::Identity(int dim)
{
	if (dim < 1)
		return false;

	this->dimX = dim;
	this->dimY = dim;

	Zero(dim, dim);
	for (int i = 0; i < dim; i++)
		this->mat[i][i] = 1;

	return true;
}

void Matrix::Show()
{
	for (int i = 0; i < dimX; i++)
	{
		for (int j = 0; j < dimY; j++)
			std::cout << this->mat[i][j] << ' ';

		std::cout << '\n';
	}
}

double Matrix::GetMean()
{
	double result = 0.0;
	for (int i = 0; i < dimX; i++)
	{
		for (int j = 0; j < dimY; j++)
			result += this->mat[i][j];
	}

	return result / (dimX + dimY);
}

Matrix Matrix::operator+(Matrix pMat)
{
	if (this->dimX != pMat.dimX || this->dimY != pMat.dimY)
	{
		std::cerr << "diff dimX or dimY" << std::endl;
		exit(-1);
	}

	Matrix res;
	res.Zero(this->dimX, this->dimY);

	for (int i = 0; i < this->dimX; i++)
	{
		for (int j = 0; j < this->dimY; j++)
			res.mat[i][j] = this->mat[i][j] + pMat.mat[i][j];
	}

	return res;
}

Matrix Matrix::operator-(Matrix pMat)
{
	if (this->dimX != pMat.dimX || this->dimY != pMat.dimY)
	{
		std::cerr << "diff dimX or dimY" << std::endl;
		exit(-1);
	}

	Matrix res;
	res.Zero(this->dimX, this->dimY);

	for (int i = 0; i < this->dimX; i++)
	{
		for (int j = 0; j < this->dimY; j++)
			res.mat[i][j] = this->mat[i][j] - pMat.mat[i][j];
	}

	return res;
}

Matrix Matrix::operator*(Matrix pMat)
{
	if(this->dimY != pMat.dimX)
	{
		std::cerr << "diff dimY and dimX" << std::endl;
		exit(-1);
	}

	Matrix res;
	res.Zero(this->dimX, pMat.dimY);

	for (int i = 0; i < this->dimX; i++)
	{
		for (int j = 0; j < pMat.dimY; j++)
		{
			for (int k = 0; k < this->dimY; k++)
				res.mat[i][j] += this->mat[i][k] * pMat.mat[k][j];
		}
	}

	return res;
}

Matrix Matrix::Transpose()
{
	Matrix res;
	res.Zero(this->dimY, this->dimX);

	for (int i = 0; i < this->dimX; i++)
	{
		for (int j = 0; j < this->dimY; j++)
			res.mat[j][i] = this->mat[i][j];
	}

	return res;
}

Matrix Matrix::Cholesky(double ztol = 1.0e-5)
{
	Matrix res;
	res.Zero(this->dimX, this->dimX);

	for (int i = 0; i < this->dimX; i++)
	{
		std::vector<double> tmp;
		double S = 0.0;
		double d = 0.0;
		for (int k = 0; k < i; k++)
			tmp.push_back(res.mat[k][i] * res.mat[k][i]);
		
		S = std::accumulate(tmp.begin(), tmp.end(), 0.0);
		d = this->mat[i][i] - S;

		if (std::abs(d) < ztol)
			res.mat[i][i] = 0.0;
		else
		{
			if (d < 0.0)
			{
				std::cerr << "less than 0.0" << std::endl;
				exit(-1);
			}

			res.mat[i][i] = std::sqrt(d);
		}

		for (int j = i + 1; j < this->dimX; j++)
		{
			std::vector<double> tmp;
			for (int k = 0; k < this->dimX; k++)
				tmp.push_back(res.mat[k][i] * res.mat[k][j]);

			S = std::accumulate(tmp.begin(), tmp.end(), 0.0);

			if (std::abs(S) < ztol)
				S = 0.0;

			try {
				res.mat[i][j] = (this->mat[i][j] - S) / res.mat[i][i];
			}
			catch (std::exception& e) {
				std::cerr << "Exception: " << e.what() << std::endl;
				exit(-1);
			}
		}
	}

	return res;
}

Matrix Matrix::CholeskyInverse()
{
	Matrix res;
	res.Zero(this->dimX, this->dimX);

	double S;
	double tjj;
	for (int j = this->dimX - 1; j >= 0; j--)
	{
		tjj = this->mat[j][j];

		std::vector<double> tmp;
		for (int k = j + 1; k < this->dimX; k++)
			tmp.push_back(this->mat[j][k] * res.mat[j][k]);

		S = std::accumulate(tmp.begin(), tmp.end(), 0.0);
		res.mat[j][j] = 1.0 / (tjj * tjj) - S / tjj;

		for (int i = j - 1; i >= 0; i--)
		{
			std::vector<double> tmp;
			for (int k = i + 1; k < this->dimX; k++)
				tmp.push_back(this->mat[i][k] * res.mat[k][j]);

			res.mat[j][i] = res.mat[i][j] = -std::accumulate(tmp.begin(), tmp.end(), 0.0) / this->mat[i][i];
		}
	}

	return res;
}

Matrix Matrix::Inverse()
{
	Matrix aux = Cholesky();
	Matrix res = aux.CholeskyInverse();

	return res;
}