#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

namespace {
	std::vector<double> operator+(const std::vector<double>& lhs, const std::vector<double>& rhs)
	{
		std::vector<double> result;
		for (int i = 0; i < (int)lhs.size(); i++)
		{
			result.push_back(lhs.at(i) + rhs.at(i));
		}
		return result;
	}

	std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& lhs, const std::vector<std::vector<double>>& rhs)
	{
		std::vector<std::vector<double>> result{ lhs };
		for (int i = 0; i < (int)lhs.size(); i++)
			for (int j = 0; j < (int)lhs.size(); j++) {
				result[i][j] += rhs[i][j];  // lhs already added
			}
		return result;
	}

	std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& lhs, const std::vector<std::vector<double>>& rhs)
	{
		std::vector<std::vector<double>> result{ lhs };
		for (int i = 0; i < (int)lhs.size(); i++)
			for (int j = 0; j < (int)lhs.size(); j++) {
				result[i][j] -= rhs[i][j];  // lhs already added
			}
		return result;
	}


	std::vector<std::vector<double>> operator*(const double& n, const std::vector<std::vector<double>>& rhs)
	{
		std::vector<std::vector<double>> result{ rhs };
		for (int i = 0; i < (int)rhs.size(); i++)
			for (int j = 0; j < (int)rhs.size(); j++) {
				result[i][j] *= n;  // rhs already added 
			}
		return result;
	}

	std::vector<double> matmul3(std::vector<std::vector<double>>& M, std::vector<double>& vec) {

		std::vector<double> results{ vec };
		results[0] = (M[0][0] * vec[0] + M[0][1] * vec[1] + M[0][2] * vec[2]);
		results[1] = (M[1][0] * vec[0] + M[1][1] * vec[1] + M[1][2] * vec[2]);
		results[2] = (M[2][0] * vec[0] + M[2][1] * vec[1] + M[2][2] * vec[2]);
		return results;
	}

	std::vector<std::vector<double>> mattranspose3(const std::vector<std::vector<double>>& M) {

		std::vector<std::vector<double>> results{ {M[0][0], M[1][0], M[2][0]},
													{M[0][1], M[1][1], M[2][1]},
													{M[0][2], M[1][2], M[2][2]} };
		return results;
	}

	std::vector<std::vector<double>> matmul3(const std::vector<std::vector<double>>& lhs, const std::vector<std::vector<double>>& rhs) {

		std::vector<std::vector<double>> results{ {0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0} };
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				for (int k = 0; k < 3; ++k) {
					results[i][j] += lhs[i][k] * rhs[k][j];
				}
			}
		}
		return results;
	}

	std::vector<std::vector<double>> matinverse3(const std::vector<std::vector<double>>& M) {
		double determinant = 0.0;
		std::vector<std::vector<double>> results{ {0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0} };
		for (int i = 0; i < 3; i++)
			determinant = determinant + (M[0][i] * (M[1][(i + 1) % 3] * M[2][(i + 2) % 3] - M[1][(i + 2) % 3] * M[2][(i + 1) % 3]));
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				results[i][j] = ((M[(j + 1) % 3][(i + 1) % 3] * M[(j + 2) % 3][(i + 2) % 3]) - (M[(j + 1) % 3][(i + 2) % 3] * M[(j + 2) % 3][(i + 1) % 3])) / determinant;
			}
		}
		return results;
	}
}

class AttitudeEstimation {

	const double N = 10.0; // propagation steps
	double t_sample = 0.0;
	double velocity = 0.0; // velocity of the system
	const double g = 9.81;

	// init matrices
	std::vector<std::vector<double>> P{
		{1e-1, 0.0, 0.0},
		{0.0, 1e-1, 0.0},
		{0.0, 0.0, 0.0}
	};
	std::vector<std::vector<double>> Q{ P };
	std::vector<std::vector<double>> R{ P };
	std::vector<std::vector<double>> A{ // f-jacobian
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0}
	};
	std::vector<std::vector<double>> L{ A }; // kalman gain
	std::vector<std::vector<double>> C{ A }; // y-jacobian
	const std::vector<std::vector<double>> I{
		{1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0}
	};

	std::vector<double> x_hat{ 0.0, 0.0, 0.0 }; // {phi_hat, theta_hat, psi_hat}
	std::vector<double> y{ 0.0, 0.0, 0.0 }; // output // accel

	double p = 0.0, q = 0.0, r = 0.0; // gyro inputs
	double phi = x_hat[0], theta = x_hat[1], psi = x_hat[2]; // estimated states


public:
	AttitudeEstimation(const double& sigma_gyro, const double& sigma_accel, const double& dt, const double& velocity = 0.0) {
		this->velocity = velocity;
		this->t_sample = dt / N;
		for (int i = 0; i < 3; i++) {
			Q[i][i] = sigma_gyro;
			R[i][i] = sigma_accel;
		}
	}

	std::vector<double> update(const double& gyro_x, const double& gyro_y, const double& gyro_z, const double& accel_x, const double& accel_y, const double& accel_z) {

		p = gyro_x, q = gyro_y, r = gyro_z;

		// propogation step
		for (int n = 0; n < N; n++) {
			x_hat[0] += t_sample * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta));
			x_hat[1] += t_sample * (q * cos(phi) - r * sin(phi));
			x_hat[2] += t_sample * (q * sin(phi) / cos(theta) + r * cos(phi) / cos(theta));

			phi = x_hat[0], theta = x_hat[1], psi = x_hat[2];

			// jacobian  -- most of them are 0
			A[0][0] = q * cos(phi) * tan(theta) - r * sin(phi) * tan(theta);
			A[0][1] = (q * sin(phi) - r * cos(phi)) / (cos(theta) * cos(theta));
			A[1][0] = -q * sin(phi) - r * cos(phi);

			P = P + t_sample * (matmul3(A, P) + matmul3(P, mattranspose3(A)) + Q);
		}

		// measurement step
		// y_ = y - h = accel input - accel equ
		y[0] = accel_x - (q * velocity * sin(theta) + g * sin(theta));
		y[1] = accel_y - (r * velocity * cos(theta) - p * velocity * sin(theta) - g * cos(theta) * sin(phi));
		y[2] = accel_z - (-q * velocity * cos(theta) - g * cos(theta) * cos(phi));

		// jacobian  -- most of them are 0 
		C[0][1] = q * velocity * cos(theta) + g * cos(theta);
		C[1][0] = -g * cos(phi) * cos(theta);
		C[1][1] = -r * velocity * cos(theta) - p * velocity * cos(theta) + g * sin(theta) * sin(phi);
		C[2][0] = g * cos(theta) * sin(phi);
		C[2][1] = q * velocity * sin(theta) + g * sin(theta) * cos(phi);

		// kalman gain
		L = matmul3(matmul3(P, mattranspose3(C)), matinverse3(R + matmul3(matmul3(C, P), mattranspose3(C))));
		x_hat = x_hat + matmul3(L, y);
		P = matmul3((I - matmul3(L, C)), P);

		return x_hat; // current estimation
	}
};
