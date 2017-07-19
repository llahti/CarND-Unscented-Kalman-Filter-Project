#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
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

  // Initialize radar measurement covariance matrix
  S_rad_ = MatrixXd(n_z_rad_, n_z_rad_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // Initialize sigmapoints matrix
  Xsig_ = MatrixXd(n_sp_x_, n_x_);
  Xsig_.fill(0.0);

  // Initialize predicted sigma points with zero
  Xsig_pred_ = MatrixXd(n_x_, n_sp_xaug_);
  Xsig_pred_.fill(0.0);

  // Initialize augmented sigma points with zero
  Xsig_aug_ = MatrixXd(n_aug_, n_sp_xaug_);
  Xsig_aug_.fill(0.0);
}

UKF::~UKF() {}

/**
 * @brief UKF::FirstUpdate handles the very first update cycle and initializes state vector.
 * @param measurement_pack
 */
void UKF::FirstUpdate(MeasurementPackage measurement_pack) {

   // Initialize state vector
   VectorXd x_tmp = VectorXd(n_x_);
   //x_tmp.fill(0.0);

   if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
     /**
     Convert radar from polar to cartesian coordinates and initialize state.
     */
     VectorXd polar = measurement_pack.raw_measurements_;
     x_tmp =  tools.Polar2Cartesian(polar);
     // Radar do not have enough information to predict speed, yaw and yaw rate so initialize them as zero
     x_tmp(2) = 0; x_tmp(3) = 0, x_tmp(4) = 0;
   }
   else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
     /**
     Initialize state.
     */
     //set the state with the initial location and zero velocity, and init speed, yaw and yaw rate as zero
     x_tmp << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0, 0.0;
   }
   previous_timestamp_ = measurement_pack.timestamp_;

   x_ = x_tmp;

   // done initializing, no need to predict or update
   is_initialized_ = true;
   return;
}

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
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  /*
  if (!is_initialized_) {
      // Init first measurement
      return FirstUpdate(meas_package);
    }
  else {
  // Predict
  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	// dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
  // Predict k+1 state
  Prediction(dt);

  // Update
  return;
  }
  */

  /* *******************
   * * TESTING METHODS *
   * *******************/
  MatrixXd Xsig;
  if (!is_initialized_) {
    Xsig_pred_ <<
             5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
               1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
              2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
             0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
              0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

    x_ <<
       5.93637,
       1.49035,
       2.20528,
      0.536853,
      0.353577;

    P_ <<
      0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
      -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
      0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
     -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
     -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

    //create example matrix with sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_rad_, n_sp_xaug_);
    Zsig <<
        6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
       0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
        2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;

    VectorXd z_pred = VectorXd(n_z_rad_);
    z_pred <<
        6.12155,
       0.245993,
        2.10313;

    S_rad_ <<
        0.0946171, -0.000139448,   0.00407016,
     -0.000139448,  0.000617548, -0.000770652,
       0.00407016, -0.000770652,    0.0180917;
    //create example vector for incoming radar measurement
    VectorXd z = VectorXd(n_z_rad_);
    z <<
        5.9214,
        0.2187,
        2.0062;


    //radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    //radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.0175;

    //radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.1;

    std_a_ = 0.2;
    std_yawdd_ = 0.2;
    //lambda_ = 3 - n_aug_;
    double delta_t = 0.1; //time diff in sec
    //cout << Xsig_pred_;
    UpdateRadarState(&Zsig, &z_pred, &z);
    is_initialized_ = true;
    //print result
    std::cout << "x_: " << std::endl << x_ << std::endl;
    std::cout << "P_: " << std::endl << P_ << std::endl;
    }

  return;

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

  // 1. Generate sigma points
  GenerateSigmaPoints(&Xsig_);
  //std::cout << "Sigma Points\n" << Xsig_ << std::endl << std::endl;
  // 2. Predict sigma points
  AugmentedSigmaPoints(&Xsig_aug_);
  //std::cout << "Augmented Sigma Points\n" << Xsig_aug_ << std::endl << std::endl;
  SigmaPointPrediction(&Xsig_pred_, &Xsig_aug_, delta_t);
  //std::cout << "Predicted Sigma Points\n" << Xsig_pred_ << std::endl << std::endl;
  // 3. Predict mean and covariance matrix
  PredictMeanAndCovariance(&Xsig_pred_);
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
  VectorXd z_pred;
  MatrixXd Zsig;
  PredictRadarMeasurement(&z_pred, &Zsig);
  UpdateRadarState(&Zsig, &z_pred, &meas_package.raw_measurements_);
}


/**
 * This method generates sigma-points and calculations depends on instance variables P_, n_x_, n_sp_x_ and x_.
 * @brief UKF::GenerateSigmaPoints
 * @param Xsig_out
 */
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, n_sp_x_);

  //set first column of sigma point matrix
  Xsig.col(0) = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }

  //write result
  *Xsig_out = Xsig;
}

/**
 * @brief UKF::AugmentedSigmaPoints
 * @param Xsig_out Augmented sigma points
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //x_aug.fill(0.0);
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sp_xaug_);
  Xsig_aug.fill(0.0);



  //create augmented covariance matrix by inserting P_ matrix to top-left corner
  P_aug.topLeftCorner(5,5) = P_;
  // Insert process noise to bottom right corner
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;  // Mean state
  for (int i = 0; i<n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_aug_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_aug_+n_aug_) * L.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug;
}


/**
 * @brief UKF::SigmaPointPrediction
 * @param Xsig_out Predicted sigma points out
 * @param Xsig_aug Augmented sigmapoints in
 * @param delta_t time difference between this and previus state
 */
void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, MatrixXd* Xsig_aug, double delta_t) {
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, n_sp_xaug_);

  //predict sigma points
  for(int i=0; i < n_sp_xaug_; i++){
    // Variables to store (previous state) xk values for easier access to them
    double p_x = (*Xsig_aug)(0, i);
    double p_y = (*Xsig_aug)(1, i);
    double v = (*Xsig_aug)(2, i);
    double yaw = (*Xsig_aug)(3, i);
    double yawd = (*Xsig_aug)(4, i);
    double nu_a = (*Xsig_aug)(5, i);
    double nu_yawdd = (*Xsig_aug)(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero when angle change is zero or very close to it
    if (fabs(yawd) > 0.0001){
      px_p = p_x + (v / yawd) * ( sin(yaw + (yawd * delta_t)) - sin(yaw));  // predict px_k+1
      py_p = p_y + (v / yawd) * ( + cos(yaw) - cos(yaw + (yawd * delta_t)));
    }
    else {
      // Predict px
      p_x = p_x + v * cos(yaw) * delta_t;
      p_y = p_y + v * sin(yaw) * delta_t;
    }

    // Predict remaining elements
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // Add noise
    px_p += 0.5 * pow(delta_t, 2) * cos(yaw) * nu_a;  // Add noise to px_k+1
    py_p += 0.5 * pow(delta_t, 2) * sin(yaw) * nu_a;
    v_p  += delta_t * nu_a;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma points into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }
  //write result
  *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(MatrixXd* Xsig_pred) {
  //create vector for weights
  VectorXd weights = VectorXd(n_sp_xaug_);

  //set weights
  // TODO: Move to instance constructor
  weights(0) = float(lambda_aug_) / float((lambda_aug_ + n_aug_));
  for (int i=1; i < n_sp_xaug_; i++){
    float w = 1. / (2 * (lambda_aug_ + n_aug_));
    weights(i) = w;
  }
  //predict state mean
  x_.fill(0.0);
  for (int i=0; i < n_sp_xaug_; i++){
    x_ += (*Xsig_pred).col(i) * weights(i);
  }
  //predict state covariance matrix
  P_.fill(0.0);
  for (int i=0; i < n_sp_xaug_; i++){
      // Calculate difference
      VectorXd x_diff = (*Xsig_pred).col(i) - x_;
      // Normalize angles. TODO: Refactor to separate function
      while (x_diff(3) > M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;
      P_ = P_ + weights(i) * x_diff * x_diff.transpose();
  }

}

/**
 * @brief UKF::PredictRadarMeasurement
 * @param z_out
 * @param Zsig_out
 *
 */
void UKF::PredictRadarMeasurement(VectorXd* z_pred_out, MatrixXd* Zsig_out) {

  //set vector for weights
  VectorXd weights = VectorXd(n_sp_xaug_);
   //double weight_0 = lambda_aug_/(lambda_aug_+n_aug);
  weights(0) = lambda_aug_/(lambda_aug_+n_aug_);
  double f_ = n_aug_+lambda_aug_;
  for (int i=1; i<n_sp_xaug_; i++) {
    double weight = 0.5/f_;
    weights(i) = weight;
  }

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
    z_pred += Zsig.col(i) * weights(i);
  }

  // R-matrix
  MatrixXd R = MatrixXd(3, 3);
  R << pow(std_radr_, 2), 0, 0,
       0, pow(std_radphi_, 2), 0,
       0, 0, pow(std_radrd_, 2);

  //calculate measurement covariance matrix S
  S_rad_.fill(0.0);
  for (int i=0; i < n_sp_xaug_; i++){
      // Calculate difference
      VectorXd z_diff = Zsig.col(i) - z_pred;
      // Normalize angles
      while (z_diff(1) > M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;
      S_rad_ = S_rad_ + weights(i) * z_diff * z_diff.transpose();
  }
  S_rad_ += R;

  //write result
  *z_pred_out = z_pred;
  *Zsig_out = Zsig;
}


void UKF::UpdateRadarState(MatrixXd* Zsig, VectorXd* z_pred, VectorXd* z)
{
  //set vector for weights
  VectorXd weights = VectorXd(n_sp_xaug_);
   double weight_0 = lambda_aug_/(lambda_aug_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<n_sp_xaug_; i++) {
    double weight = 0.5/(n_aug_+lambda_aug_);
    weights(i) = weight;
  }

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_rad_);

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

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();

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

  //write result
  //*x_out = x_;
  //*P_out = P_;
}

