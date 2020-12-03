#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

// using Eigen::MatrixXd;
// using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Dimensions
  n_x_ = 5;

  // initial state vector
  x_ = Eigen::VectorXd(n_x_);

  // initial covariance matrix
  P_ = Eigen::MatrixXd(n_x_, n_x_);

  // Augmented dimension
  n_aug_ = n_x_ + 2;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  lambda_ = 3 - n_aug_;

  size_ = 2 * n_aug_ + 1;

  weights_ = Eigen::VectorXd(size_);

  // AJUST VALUE ON weights_
  double weight_0 = lambda_ / (lambda_+n_aug_);
  double weight = 0.5/(lambda_+n_aug_);
  weights_(0) = weight_0;

  for (int i = 1; i < size_; ++i) {
      weights_(i) = weight;
  }

  // GENERATE SOME VETOR AND MATRIX's
  x_aug = Eigen::VectorXd(n_aug_);
  P_aug = Eigen::MatrixXd(n_aug_, n_aug_);
  Xsig_aug = Eigen::MatrixXd(n_aug_, size_);
  Xsig_pred_ = Eigen::MatrixXd(n_x_, size_);
  Q_ = Eigen::MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
  Q_ << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;
  mea_sig_ = Eigen::MatrixXd(n_aug_ - n_x_, size_);
  mean_mea_ = Eigen::VectorXd(n_aug_ - n_x_);
  mea_cov_ = Eigen::MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
  
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
 
  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << meas_package.raw_measurements_[0], 
            meas_package.raw_measurements_[1], 
                                            0, 
                                            0, 
                                            0;
    } else {
      // COORDINATE TRANSFORMATION POLAR TO RECT
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double x = rho * std::cos(phi);
      double y = rho * std::sin(phi);
      double vx = rho_dot * std::cos(phi);
      double vy = rho_dot * std::sin(phi);
      double v = std::sqrt(vx * vx + vy * vy);
      x_ << x, y, v, 0, 0;
    }
    is_initialized_ = true;
    return;
  }
  //std::cout << meas_package.raw_measurements_ << std::endl;

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  UKF::Prediction(dt);

  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UKF::UpdateLidar(meas_package);
  }

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UKF::UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // PREDICTION
  x_aug.fill(0.0);
  x_aug.head(5) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  // PREDICTION COVARIANCE MATRIX
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q_;

  // SQUARE ROOT MATRIX
  SRM_ = P_aug.llt().matrixL();

  // SIGMA POINTS
  Xsig_aug.fill(0.0);
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + std::sqrt(lambda_ + n_aug_) * SRM_.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - std::sqrt(lambda_+n_aug_) * SRM_.col(i);
  }

  // PREDICTION OF SIGMA POINTS
  for (int i = 0; i < size_; i++)
  {
    double _x = Xsig_aug(0, i);
    double _y = Xsig_aug(1, i);
    double _v = Xsig_aug(2, i);
    double _yaw = Xsig_aug(3, i);
    double _yawd = Xsig_aug(4, i);
    double _noise = Xsig_aug(5, i);
    double _noise_yaw = Xsig_aug(6, i);

    if (std::fabs(_yawd) > 0.001)
    {
      Xsig_pred_(0, i) = _x + _v / _yawd * (std::sin(_yaw + _yawd * delta_t) - std::sin(_yaw));
      Xsig_pred_(1, i) = _y + _v / _yawd * (-1) * (std::cos(_yaw + _yawd * delta_t) - std::cos(_yaw));
    } else {
      Xsig_pred_(0, i) = _x + _v * delta_t * std::cos(_yaw);
      Xsig_pred_(1, i) = _y + _v * delta_t * std::sin(_yaw);
    }

    // ADDING SOME NOISE
    Xsig_pred_(0, i) = Xsig_pred_(0, i) + 0.5 * _noise * delta_t * delta_t * std::cos(_yaw);
    Xsig_pred_(1, i) = Xsig_pred_(1, i) + 0.5 * _noise * delta_t * delta_t * std::sin(_yaw); 
    Xsig_pred_(2, i) = _v + _noise * delta_t;
    Xsig_pred_(3, i) = _yaw + _yawd * delta_t + 0.5 * _noise_yaw * delta_t * delta_t;
    Xsig_pred_(4, i) = _yawd + _noise_yaw * delta_t;

  }

  // PREDICTION OF THE STATE MEAN
  x_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
  //PREDICTION OF STATE COVARIANCE
  P_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //NORMALIZATION
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }


}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  mea_ = meas_package.raw_measurements_;
  for (int i = 0; i < size_; i++)
  {
    mea_sig_(0, i) = Xsig_pred_(0, i);
    mea_sig_(1, i) = Xsig_pred_(1, i);
  }
  // Prediction of the mean measurement
  mean_mea_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    mean_mea_ += weights_(i) * mea_sig_.col(i);
  }
  // Calculate covariance
  mea_cov_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    Eigen::VectorXd mea_diff_ = mea_sig_.col(i) - mean_mea_;
    mea_cov_ += weights_(i) * mea_diff_ * mea_diff_.transpose();
  }

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}