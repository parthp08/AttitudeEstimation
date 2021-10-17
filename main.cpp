#include <iostream>
#include <vector>
#include <fstream>
#include "AttitudeEstimation.hpp"
#include "ReadCSV.hpp"

constexpr double RAD2DEG = 57.2958;

int main()
{
    // reading imu_data collected from phone
    auto imu_data = read_csv("imu_data.csv");
    std::vector<double> t = std::get<1>(imu_data[0]);
    std::vector<double> accel_x = std::get<1>(imu_data[1]);
    std::vector<double> accel_y = std::get<1>(imu_data[2]);
    std::vector<double> accel_z = std::get<1>(imu_data[3]);
    std::vector<double> gyro_x = std::get<1>(imu_data[4]);
    std::vector<double> gyro_y = std::get<1>(imu_data[5]);
    std::vector<double> gyro_z = std::get<1>(imu_data[6]);
    
    // storing output for plot
    std::vector<double> euler_angles{ 0.0, 0.0, 0.0 };
    std::vector<double> phi{ t };
    std::vector<double> theta{ t };
    std::vector<double> psi{ t };

    // gyro and accel sensor varinace // from datasheets
    const double sigma_gyro = 1e-20;
    const double sigma_accel = 0.0025 * 0.0025;

    // sampling time of sensors
    const double dt = t[1] - t[0];
    
    AttitudeEstimation attitude(sigma_gyro, sigma_accel, dt); // init

    for (int i = 0; i < (int)t.size(); i++) {
        euler_angles = attitude.update(gyro_x[i], gyro_y[i], gyro_z[i], accel_x[i], accel_y[i], accel_z[i]);
        phi[i] = euler_angles[0]*RAD2DEG;
        theta[i] = euler_angles[1]*RAD2DEG;
        psi[i] = euler_angles[2]*RAD2DEG;
    }

    std::ofstream output("results.csv");
    output << "t" << "," << "phi" << "," << "theta" << "," << "psi" << std::endl;
    for (int i = 0; i < t.size(); i++) {
        output << t[i] << "," << phi[i] << "," << theta[i] << "," << psi[i] << std::endl;
    }
    output.close();

    return 0;
}
