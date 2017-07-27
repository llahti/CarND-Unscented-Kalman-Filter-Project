#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;



UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Initial state is "not initialized". UKF is initialized properly during the very first measurement update.
  is_initialized_ = false;

  //
  previous_timestamp_ = 0;

  // Definitions for matrix sizes
  n_x_ = 5;  // Length of state matrix
  n_aug_ = 7;  // Length of augmented matrix
  n_sp_xaug_ =  2 * n_aug_ + 1;  // Number of sigma point vectors for augmented state matrix
  n_sp_x_ = 2 * n_x_ + 1;  // Number of sigma point vectors for state matrix
  lambda_ = double(3 - n_x_); // define spreading parameter lambda
  lambda_aug_ = 3 - n_aug_; // define spreading parameter lambda for augmented sigma-points
  n_z_rad_ = 3;  // Radar returns 3 measurement values: distance, speed, angle
  n_z_lidar_ = 2;  // Lidar returns 2 measurement values: px and py


  // initial state vector
  // px, py, velocity, yaw
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1., 0., 0., 0., 0.,
        0., 1., 0., 0., 0.,
        0., 0., 1., 0., 0.,
        0., 0., 0., 1., 0.,
        0., 0., 0., 0., 1.;

  // Initialize covariance matrix for augmented sigmapoints
  P_aug_ = MatrixXd(n_aug_, n_aug_);
  P_aug_.fill(0.0);

  // Initialize radar measurement covariance matrix
  S_rad_ = MatrixXd(n_z_rad_, n_z_rad_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 1;

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


  // Kalman Filter Matrices for Lidar update
  H_lidar_ = MatrixXd(2, 5);
  H_lidar_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  // Initialize lidar measurement covariance matrix
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << pow(std_laspx_, 2), 0,
              0, pow(std_laspy_, 2);

  // Initialize radar measurement covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << pow(std_radr_, 2), 0, 0,
              0, pow(std_radphi_, 2), 0,
              0, 0, pow(std_radrd_, 2);

  // Initialize sigmapoints matrix
  Xsig_ = MatrixXd(n_x_, n_sp_x_);
  Xsig_.fill(0.0);

  // Initialize predicted sigma points with zero
  Xsig_pred_ = MatrixXd(n_x_, n_sp_xaug_);
  Xsig_pred_.fill(0.0);

  // Initialize augmented sigma points with zero
  Xsig_aug_ = MatrixXd(n_aug_, n_sp_xaug_);
  Xsig_aug_.fill(0.0);

  // Initialize weights
  InitWeights();

  // Initialize NIS
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
}

UKF::~UKF() {}

/**
 * Initializes Unscented Kalman filter
 */
void UKF::InitWeights()
{
  weights_ = VectorXd(n_sp_xaug_);
  weights_(0) = lambda_aug_/(lambda_aug_+n_aug_);
  double f_ = n_aug_+lambda_aug_;
  for (int i=1; i<n_sp_xaug_; i++) {
    //double weight = 0.5/f_;
    weights_(i) = 0.5/f_;
  }
}

/**
 * @brief UKF::FirstUpdate handles the very first update cycle and initializes state vector.
 * @param measurement_pack
 */
void UKF::FirstUpdate(MeasurementPackage measurement_pack) {
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
    */
    double rho = measurement_pack.raw_measurements_(0);  // range
    double phi = measurement_pack.raw_measurements_(1);  // angle between rho and x
    //double rho_dot = measurement_pack.raw_measurements_(2);  // change of p (change rate)

    // Radar do not have enough information to predict
    // speed, yaw and yaw rate so initialize them as zero
    x_ << rho * cos(phi), rho * sin(phi), 0.0, 0.0, 0.0;

    // Initialize P_ matrix with measurement type specific values
    P_ << 0.5, 0., 0., 0., 0.,
          0., 0.5, 0., 0., 0.,
          0., 0., 1., 0., 0.,
          0., 0., 0., 1., 0.,
          0., 0., 0., 0., 1.;
   }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    /**
    Initialize state.
    */
    //set the state with the initial location and zero velocity, and init speed, yaw and yaw rate as zero
    x_ << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1],
          0.0,
          0.0,
          0.0;
    // Initialize P_ matrix with measurement type specific values
    P_ << 0.5, 0., 0., 0., 0.,
          0., 0.5, 0., 0., 0.,
          0., 0., 1., 0., 0.,
          0., 0., 0., 1., 0.,
          0., 0., 0., 0., 1.;
  }
  previous_timestamp_ = measurement_pack.timestamp_;

  // Initialization done, save first measurement as a current state.
  //x_ = x_tmp;
  is_initialized_ = true;
  return;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /******************
   *  Initialization
   ******************/
  //std::cout << "Meas: " << meas_package.raw_measurements_ << std::endl;
  if (!is_initialized_) {
    // Init first measurement
    FirstUpdate(meas_package);
    return;
  }

  /***********
   * Predict
   ***********/
  //compute the time in seconds elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	// dt - expressed in seconds
  // Predict k+1 state
  Prediction(dt);

  /**********
   * Update
   **********/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true) {
    // Radar updates
    // update timestamp only when measurement is processed
    previous_timestamp_ = meas_package.timestamp_;
    UpdateRadar(meas_package);
    }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true) {
    // Lidar update
    // update timestamp only when measurement is processed
    previous_timestamp_ = meas_package.timestamp_;
    UpdateLidar(meas_package);
  }
  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // 1. Generate sigma points
  GenerateSigmaPoints();
  // 2. Predict sigma points
  AugmentedSigmaPoints();
  SigmaPointPrediction(delta_t);
  // 3. Predict mean and covariance matrix
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // Predict LIDAR measurement by using sigma-points
  VectorXd z_pred;
  PredictLidarMeasurement(&z_pred);
  //VectorXd y = z - z_pred;
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  MatrixXd Ht = H_lidar_.transpose();
  MatrixXd S = H_lidar_ * P_ * Ht + R_lidar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  //long x_size = x_.size();
  MatrixXd I_ = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I_ - K * H_lidar_) * P_;

  // Calculate NIS
  NIS_laser_ = y.transpose() * Si * y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  VectorXd z_pred;
  MatrixXd Zsig;
  PredictRadarMeasurement(&z_pred, &Zsig);
  UpdateRadarState(&Zsig, &z_pred, &meas_package.raw_measurements_);

  // Calculate NIS
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  NIS_radar_ = y.transpose() * S_rad_.inverse() * y;
}


/**
 * This method generates sigma-points and calculations depends on instance variables P_, n_x_, n_sp_x_ and x_.
 * @brief UKF::GenerateSigmaPoints
 * @param Xsig_out
 */
void UKF::GenerateSigmaPoints() {
  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  Xsig_.fill(0.0);
  //set first column of sigma point matrix
  Xsig_.col(0) = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig_.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
    Xsig_.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }
}

/**
 * @brief UKF::AugmentedSigmaPoints
 */
void UKF::AugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented state covariance
  P_aug_.fill(0.0);

  //create sigma point matrix
  Xsig_aug_.fill(0.0);

  //create augmented covariance matrix by inserting P_ matrix to top-left corner
  P_aug_.topLeftCorner(5,5) = P_;
  // Insert process noise to bottom right corner
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug;  // Mean state
  for (int i = 0; i<n_aug_; i++)
  {
    Xsig_aug_.col(i+1)        = x_aug + sqrt(lambda_aug_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_aug_+n_aug_) * L.col(i);
  }
}


/**
 * @brief UKF::SigmaPointPrediction
 * @param delta_t time difference between this and previus state
 */
void UKF::SigmaPointPrediction(double delta_t) {
  Xsig_pred_.fill(0.0);

  //predict sigma points
  for(int i=0; i < n_sp_xaug_; i++){
    // Variables to store (previous state) xk values for easier access to them
    double p_x = Xsig_aug_(0, i);
    double p_y = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero when angle change is zero or very close to it
    if (fabs(yawd) > 0.0001){
      px_p = p_x + (v / yawd) * ( sin(yaw + (yawd * delta_t)) - sin(yaw));  // predict px_k+1
      py_p = p_y + (v / yawd) * ( + cos(yaw) - cos(yaw + (yawd * delta_t)));
    }
    else {
      // Predict px
      px_p = p_x + (v * delta_t * cos(yaw));
      py_p = p_y + (v * delta_t * sin(yaw));
    }

    // Predict remaining elements
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // Add noise
    px_p = px_p + 0.5 * pow(delta_t, 2) * cos(yaw) * nu_a;  // Add noise to px_k+1
    py_p = py_p + 0.5 * pow(delta_t, 2) * sin(yaw) * nu_a;
    v_p  += delta_t * nu_a;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma points into i-th column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance() {
  //predict state mean
  x_.fill(0.0);
  for (int i=0; i < n_sp_xaug_; i++){
    x_ += Xsig_pred_.col(i) * weights_(i);
  }
  //predict state covariance matrix
  P_.fill(0.0);
  for (int i=0; i < n_sp_xaug_; i++){
      // Calculate difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      // Normalize angles. TODO: Refactor to separate function
      //std::cout << "angle before norm:" << x_diff(3) << std::endl;
      while (x_diff(3) >  M_PI) x_diff(3) -= 2.*M_PI;
      while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
      //std::cout << "angle after norm: " << x_diff(3) << std::endl;
      P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * @brief UKF::PredictLidarMeasurement
 * @param z_pred_out
 */
void UKF::PredictLidarMeasurement(VectorXd *z_pred_out) {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_lidar_, n_sp_xaug_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_lidar_);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sp_xaug_; i++) {  //2n+1 simga points
    // measurement model
    Zsig(0,i) = Xsig_pred_(0,i);  // px
    Zsig(1,i) = Xsig_pred_(1,i);  // py
  }

  //calculate mean predicted measurement
  //predict state mean
  z_pred.fill(0.0);
  for (int i=0; i < n_sp_xaug_; i++){
    z_pred += Zsig.col(i) * weights_(i);
  }

  //write result
  *z_pred_out = z_pred;
}

/**
 * @brief UKF::PredictRadarMeasurement
 * @param z_out
 * @param Zsig_out
 *
 */
void UKF::PredictRadarMeasurement(VectorXd* z_pred_out, MatrixXd* Zsig_out) {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_rad_, n_sp_xaug_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_rad_);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sp_xaug_; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
  //calculate mean predicted measurement
  //predict state mean
  z_pred.fill(0.0);
  for (int i=0; i < n_sp_xaug_; i++){
    z_pred += Zsig.col(i) * weights_(i);
  }

  //calculate measurement covariance matrix S
  S_rad_.fill(0.0);
  for (int i=0; i < n_sp_xaug_; i++){
      // Calculate difference
      VectorXd z_diff = Zsig.col(i) - z_pred;
      // Normalize angles
      while (z_diff(1) > M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;
      S_rad_ = S_rad_ + weights_(i) * z_diff * z_diff.transpose();
  }
  S_rad_ += R_radar_;

  //write result
  *z_pred_out = z_pred;
  *Zsig_out = Zsig;
}


void UKF::UpdateRadarState(MatrixXd* Zsig, VectorXd* z_pred, VectorXd* z)
{
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_rad_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i=0; i < n_sp_xaug_; i++) {
    // Calculate difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    VectorXd z_diff = (*Zsig).col(i) - (*z_pred);
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S_rad_.inverse();

  // Residual
  VectorXd z_diff = (*z) - (*z_pred);
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_rad_*K.transpose();
}

