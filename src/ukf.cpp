#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Eigen;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.


  // Dimensions
  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2 * n_x_ + 1;
  n_aug_sig_ = 2 * n_aug_ + 1;
  n_z_laser_ = 2;
  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_radar_ = 3;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd::Zero(n_aug_sig_);

  // Predicted Sigma points
  Xsig_pred_ = MatrixXd::Zero(n_x_, n_aug_sig_);
  //create sigma point matrix
  Xsig_ = MatrixXd::Zero(n_x_, n_sig_);
  //create augmented mean vector
  x_aug_ = VectorXd::Zero(n_aug_);
  //create augmented state covariance
  P_aug_ = MatrixXd::Zero(n_aug_, n_aug_);
  //create sigma point matrix
  Xsig_aug_ = MatrixXd::Zero(n_aug_, n_aug_sig_);


  //create matrix for sigma points in measurement space
  Z_radar_sig_ = MatrixXd::Zero(n_z_radar_, n_aug_sig_);
  z_radar_pred = VectorXd::Zero(n_z_radar_);

  R_radar_ = MatrixXd::Zero(n_z_radar_, n_z_radar_);
  S_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  //create matrix for cross correlation Tc_radar_
  Tc_radar_ = MatrixXd(n_x_, n_z_radar_);
}

UKF::~UKF() = default;

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    previous_timestamp_ = meas_package.timestamp_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ << meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]),
          meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]),
          0,
          0,
          0;
      std::cout << meas_package.raw_measurements_ << std::endl << std::endl << std::endl;
      std::cout << x_ << std::endl << std::endl << std::endl << std::endl;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds

  // skip second prediction is measurements arrive at nearly the same time
  if (dt > 0.01) {
    previous_timestamp_ = meas_package.timestamp_;

    Prediction(dt);
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
   * Generate Sigma Points
   */
  //calculate square root of P
  MatrixXd SR_P = P_.llt().matrixL();
  //calculate sigma points ...
  MatrixXd A = std::sqrt(lambda_ + n_x_) * SR_P;
  Xsig_ << x_, A.colwise() + x_, (-A).colwise() + x_;

  /**
   * Augmentation
   */
  //create augmented mean state
  x_aug_ << x_, 0, 0;
  //create augmented covariance matrix
  P_aug_ << P_, MatrixXd::Zero(5, 2), MatrixXd::Zero(2, 5), (MatrixXd(2, 2) << std_a_ * std_a_, 0, 0, std_yawdd_
      * std_yawdd_).finished();
  //create square root matrix
  MatrixXd SR_P_aug = P_aug_.llt().matrixL();
  MatrixXd B = std::sqrt(lambda_ + n_aug_) * SR_P_aug;
  //create augmented sigma points
  Xsig_aug_ << x_aug_, B.colwise() + x_aug_, (-B).colwise() + x_aug_;

  /**
   * Sigma Point Prediction
   */
  for (int i = 0; i < n_aug_sig_; i++) {
    VectorXd x = Xsig_aug_.col(i);
    if (x(4) != 0) {
      Xsig_pred_.col(i) << x(0) + (x(2) / x(4)) * (sin(x(3) + x(4) * delta_t) - sin(x(3)))
          + 0.5 * (delta_t * delta_t) * cos(x(3)) * x(5),
          x(1) + (x(2) / x(4)) * (-cos(x(3) + x(4) * delta_t) + cos(x(3)))
              + 0.5 * (delta_t * delta_t) * sin(x(3)) * x(5),
          x(2) + 0 + delta_t * x(5),
          x(3) + x(4) * delta_t + 0.5 * (delta_t * delta_t) * x(6),
          x(4) + 0 + delta_t * x(6);
    } else {
      Xsig_pred_.col(i) << x(0) + x(2) * cos(x(4) * delta_t) + 0.5 * (delta_t * delta_t) * cos(x(3)) * x(5),
          x(1) + x(2) * sin(x(4) * delta_t) + 0.5 * (delta_t * delta_t) * sin(x(3)) * x(5),
          x(2) + 0 + delta_t * x(5),
          x(3) + 0 + 0.5 * (delta_t * delta_t) * x(6),
          x(4) + 0 + delta_t * x(6);
    }
  }

  /**
   * Mean and Covariance
   */

  // set weights
  weights_.setConstant(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  //predicted state mean
  x_ = Xsig_pred_ * weights_;

  //predicted state covariance matrix
  P_.fill(0.0);
  MatrixXd X_diff = Xsig_pred_.colwise() - x_;
  for (int i = 0; i < n_aug_sig_; i++) {
    //angle normalization
    while (X_diff.col(i)(3) > M_PI) X_diff.col(i)(3) -= 2. * M_PI;
    while (X_diff.col(i)(3) < -M_PI) X_diff.col(i)(3) += 2. * M_PI;
  }
  P_ = X_diff * weights_.asDiagonal() * X_diff.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  MatrixXd H_laser_ = MatrixXd(n_z_laser_, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0;
  MatrixXd R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;

  VectorXd z = VectorXd(n_z_laser_);
  z << meas_package.raw_measurements_[0],
      meas_package.raw_measurements_[1];

  VectorXd z_pred = H_laser_ * x_;
  VectorXd z_diff = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * z_diff);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;

  NIS_laser_.push_back(z_diff.transpose() * S.inverse() * z_diff);
  std::vector<float> x95;
  std::vector<float> x90;
  std::vector<float> x10;
  std::vector<float> x05;
  for (float &it : NIS_laser_) {
    if(it > 0.103){
      x95.push_back(it);
    }
    if (it > 0.211){
      x90.push_back(it);
    }
    if (it > 4.605){
      x10.push_back(it);
    }
    if (it > 5.991){
      x05.push_back(it);
    }
  }
  std::cout << "Laser NIS:" << std::endl;
  std::cout << (float)(x95.size())/NIS_laser_.size() << "\t";
  std::cout << (float)(x90.size())/NIS_laser_.size() << "\t";
  std::cout << (float)(x10.size())/NIS_laser_.size() << "\t";
  std::cout << (float)(x05.size())/NIS_laser_.size() << std::endl << std::endl << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * Predict Radar Measurement
   */
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x = Xsig_pred_.col(i);
    Z_radar_sig_.col(i) << sqrt(x(0) * x(0) + x(1) * x(1)),
        atan2(x(1), x(0)),
        (x(0) * cos(x(3)) * x(2) + x(1) * sin(x(3)) * x(2)) / sqrt(x(0) * x(0) + x(1) * x(1));
    while (Z_radar_sig_(1, i) > M_PI) Z_radar_sig_(1, i) -= 2. * M_PI;
    while (Z_radar_sig_(1, i) < -M_PI) Z_radar_sig_(1, i) += 2. * M_PI;
  }

  //mean predicted measurement
  z_radar_pred = Z_radar_sig_ * weights_;


  //add measurement noise covariance matrix
  R_radar_ << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;

  //innovation covariance matrix S_radar_
  S_radar_.fill(0.0);

  MatrixXd Z_diff = Z_radar_sig_.colwise() - z_radar_pred;
  for (int i = 0; i < n_aug_sig_; i++) {  //2n+1 simga points
    //angle normalization
    while (Z_diff.col(i)(1) > M_PI) Z_diff.col(i)(1) -= 2. * M_PI;
    while (Z_diff.col(i)(1) < -M_PI) Z_diff.col(i)(1) += 2. * M_PI;
  }
  S_radar_ = Z_diff * weights_.asDiagonal() * Z_diff.transpose() + R_radar_;


  /**
   * Update
   */
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z_radar_);
  z << meas_package.raw_measurements_[0],
      meas_package.raw_measurements_[1],
      meas_package.raw_measurements_[2];

  //calculate cross correlation matrix
  Tc_radar_.fill(0.0);

  MatrixXd X_diff = Xsig_pred_.colwise() - x_;
  // Z_diff = Z_radar_sig_.colwise()-z_radar_pred;
  for (int i = 0; i < n_aug_sig_; i++) {  //2n+1 simga points
    //angle normalization
    while (Z_diff.col(i)(1) > M_PI) Z_diff.col(i)(1) -= 2. * M_PI;
    while (Z_diff.col(i)(1) < -M_PI) Z_diff.col(i)(1) += 2. * M_PI;

    while (X_diff.col(i)(3) > M_PI) X_diff.col(i)(3) -= 2. * M_PI;
    while (X_diff.col(i)(3) < -M_PI) X_diff.col(i)(3) += 2. * M_PI;
  }
  Tc_radar_ = X_diff * weights_.asDiagonal() * Z_diff.transpose();

  //Kalman gain K;
  MatrixXd K = Tc_radar_ * S_radar_.inverse();

  //residual
  VectorXd z_diff = z - z_radar_pred;

  //angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_radar_ * K.transpose();

  NIS_radar_.push_back(z_diff.transpose() * S_radar_.inverse() * z_diff);
  std::vector<float> x95;
  std::vector<float> x90;
  std::vector<float> x10;
  std::vector<float> x05;
  for (float &it : NIS_radar_) {
    if(it > 0.352){
      x95.push_back(it);
    }
    if (it > 0.584){
      x90.push_back(it);
    }
    if (it > 6.251){
      x10.push_back(it);
    }
    if (it > 7.815){
      x05.push_back(it);
    }
  }
  std::cout << "Radar NIS:" << std::endl;
  std::cout << (float)(x95.size())/NIS_radar_.size() << "\t";
  std::cout << (float)(x90.size())/NIS_radar_.size() << "\t";
  std::cout << (float)(x10.size())/NIS_radar_.size() << "\t";
  std::cout << (float)(x05.size())/NIS_radar_.size() << std::endl << std::endl << std::endl;
}
