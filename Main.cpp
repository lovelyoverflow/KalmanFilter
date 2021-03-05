#include <iostream>
#include <cstdio>
#include <vector>
#include <initializer_list>
#include <fstream>
#include <string>
#include "Matrix.h"

std::pair<std::vector<double>, std::vector<double>> KalmanFilter(Matrix measurements, 
	Matrix x, Matrix P)
{
	std::vector<double> time;
	std::vector<double> alt;

	Matrix F(
		{
		{1., 0., 0.1, 0.},
		{0., 1., 0.,  0.1},
		{0., 0., 1.,  0.},
		{0., 0., 0.,  1.}
		}
	);
	Matrix H(
		{
		{1., 0., 0., 0.},
		{0., 1., 0., 0.}
		}
	);
	Matrix R(
		{
		{0.1, 0.},
		{0., 0.1}
		}
	);
	Matrix I(
		{
		{1., 0., 0., 0.},
		{0., 1., 0., 0.},
		{0., 0., 1., 0.},
		{0., 0., 0., 1.}
		}
	);
	Matrix u(
		{
		{0.},
		{0.},
		{0.},
		{0.}
		}
	);

	for (int i = 0; i < measurements.mat.size(); i++)
	{
		x = (F * x) + u;
		P = F * P * F.Transpose();

		Matrix z({ measurements.mat[i] });
		Matrix y(z.Transpose() - (H * x));
		Matrix S(H * P * H.Transpose() + R);
		Matrix K(P * H.Transpose() * S.Inverse());
		x = x + (K * y);
		P = (I - (K * H)) * P;

		time.push_back(x.mat[0][0]);
		alt.push_back(x.mat[1][0]);
	}

	return std::make_pair(time, alt);
}

int main()
{
	std::ifstream altDataFile;
	altDataFile.open("test2.txt");

	std::vector<std::vector<double>> measurements;
	std::vector<double> initial_xy;

	if (altDataFile.is_open())
	{
		double time = 0;
		while (!altDataFile.eof())
		{
			char str[256];
			altDataFile.getline(str, 256);
			
			if (time == 0.0)
				initial_xy = { time, atof(str) };
			else
				measurements.push_back({ time, atof(str) });
			time += 1.0;
		}
	}

	Matrix x({ { initial_xy[0]}, {initial_xy[1]}, {0.}, {0.} });
	double dt = 0.1;

	altDataFile.close();

	Matrix P(
		{
		{0., 0., 0.,    0.},
		{0., 0., 0.,    0.},
		{0., 0., 1000., 0.},
		{0., 0., 0.,    1000.}
		}
	);

	std::pair<std::vector<double>, std::vector<double>> filteredData = KalmanFilter(measurements, x, P);
	std::vector<double> timeArr(filteredData.first);
	std::vector<double> altArr(filteredData.second);

	double alt = altArr[0];
	double result = 0.0;

	for (int i = 0; i < timeArr.size(); i++)
	{
		double beforeAlt = alt;
		alt = altArr[i];

		double error = alt - beforeAlt;

		if (error > 0)
			result += error;
	}

	std::cout << result << std::endl;
	return 0;
}