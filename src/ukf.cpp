#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
#define ZERO_VAL 0.001
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
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ =  2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;
  
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
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;
  n_x_ = x_.size();
  n_aug_ = n_x_ + 2;
  lambda_ = 3-n_aug_;
  n_sig_ = 2*n_aug_ + 1;
  weights_ = VectorXd(n_sig_);
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  //initialize weights
    weights_(0) = lambda_/(lambda_+n_aug_);;
    for (int i=1; i< n_sig_; i++) {  //2n+1 weights
      double weight = 0.5/(n_aug_+lambda_);
      weights_(i) = weight;
    }

    R_Lidar_ = MatrixXd(2,2);
    R_Lidar_ <<  std_laspx_*std_laspx_, 0,
                 0, std_laspy_*std_laspy_;
    
    R_Radar_ = MatrixXd(3,3);
    R_Radar_ <<  std_radr_*std_radr_, 0, 0,
                 0, std_radphi_*std_radphi_, 0,
                 0, 0,std_radrd_*std_radrd_;
} 

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(!is_initialized_){
    time_us_ = meas_package.timestamp_;
    
    
    P_ <<  1,   0,    0,   0,   0,
           0,   1,    0,   0,   0,
           0,   0,    1,   0,   0,
           0,   0,    0,   1,   0,
           0,   0,    0,   0,   1;
    
    
    if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      //Extract values for px, py. Set velocities and rate to 0
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      //Check for values close to zero
      if(fabs(x_(0)) < ZERO_VAL and fabs(x_(1)) < ZERO_VAL){
        x_(0) = ZERO_VAL;
        x_(1) = ZERO_VAL;
      }
    }
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      //Extract values for rho, rho_dot, and phi and convert them to cartesian coords to get values for px, py, vx, vy
      const float px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      const float py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
      const float vx = meas_package.raw_measurements_[2] * cos(meas_package.raw_measurements_[1]);
      const float vy = meas_package.raw_measurements_[2] * sin(meas_package.raw_measurements_[1]);
      const float v = sqrt(vx*vx + vy*vy);
      x_ << px, py, v, 0, 0;
    }
    
  
    
    is_initialized_ = true;
    return;
  }

  const double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  if(meas_package.sensor_type_ == MeasurementPackage::LASER ){
    UpdateLidar(meas_package);
  }
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR ){
    UpdateRadar(meas_package);
  }
  
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out){
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);


  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
 
  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  //x_aug(5) = 0;
  //x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  const double lambda_aug_sqrt = sqrt(lambda_ + n_aug_);
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + lambda_aug_sqrt * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  *Xsig_out = Xsig_aug;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  //Generate augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  AugmentedSigmaPoints(&Xsig_aug);

  //Predict sigma points  
  SigmaPointPrediction(Xsig_aug, delta_t);

  //Predict mean & covariance
  x_ = Xsig_pred_ * weights_;
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    NormalizeAngle(&x_diff(3));
    

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}


void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t){
  for (int i = 0; i< n_sig_; i++)
  {
    //extract values for better readability
    const double p_x = Xsig_aug(0,i);
    const double p_y = Xsig_aug(1,i);
    const double v = Xsig_aug(2,i);
    const double yaw = Xsig_aug(3,i);
    const double yawd = Xsig_aug(4,i);
    const double nu_a = Xsig_aug(5,i);
    const double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > ZERO_VAL) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2;
  MatrixXd Zsig = Xsig_pred_.block(0,0,n_z,n_sig_);
  UpdateUKFState(meas_package,Zsig,n_z);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    int n_z = 3;
    MatrixXd Zsig = MatrixXd(n_z, n_sig_);
    //transform sigma points into measurement space
    for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

      // extract values for better readibility
      const double p_x = Xsig_pred_(0,i);
      const double p_y = Xsig_pred_(1,i);
      const double v  = Xsig_pred_(2,i);
      const double yaw = Xsig_pred_(3,i);

      const double v1 = cos(yaw)*v;
      const double v2 = sin(yaw)*v;

      // measurement model
      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
      Zsig(1,i) = atan2(p_y,p_x);                                 //phi
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    UpdateUKFState(meas_package,Zsig, n_z);

}

void UKF::UpdateUKFState(MeasurementPackage meas_package, MatrixXd Zsig, int n_z){
  //mean predicted measurement
  
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sig_; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    NormalizeAngle(&z_diff(1));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  if(meas_package.sensor_type_ == MeasurementPackage::LASER ){
    S = S + R_Lidar_;
  }
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR ){
    S = S + R_Radar_;
  }

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization if Radar
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){  
      NormalizeAngle(&z_diff(1));
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    NormalizeAngle(&x_diff(3));
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //angle normalization for Radar
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){  
    NormalizeAngle(&z_diff(1));    
  }

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}

void UKF::NormalizeAngle(double *angle){
  while (*angle > M_PI) *angle -= 2.*M_PI;
  while (*angle < -M_PI) *angle += 2.*M_PI;
}
