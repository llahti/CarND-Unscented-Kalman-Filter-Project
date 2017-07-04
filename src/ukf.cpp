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

  // Definitions for matrix sizes
  n_x_ = 5;  // Length of state matrix
  n_aug_ = 7;  // Length of augmented matrix
  n_sp_xaug_ =  2 * n_aug_ + 1;  // Number of sigma point vectors for augmented state matrix
  n_sp_x_ = 2 * n_x_ + 1;  // Number of sigma point vectors for state matrix

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1., 0., 0., 0., 0.,
        0., 1., 0., 0., 0.,
        0., 0., 1., 0., 0.,
        0., 0., 0., 1., 0.,
        0., 0., 0., 0., 1.;

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

  // Initialize predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_spoints_);
  Xsig_pred_.fill(0.0);

  Xsig_aug = MatrixXd(n_aug_, n_spoints_);
  Xsig_aug_.fill(0.0);


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
}


/**
 * @brief UKF::GenerateSigmaPoints
 * @param Xsig_out
 * @param x_in State matrix X in
 * @param P_in State covarience matrix in
 * @param lambda
 * @param n_x Number of points in state matrix
 * @param n_aug Number of sigmapoints
 */
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out, const VectorXd* x_in, const MatrixXd* P_in, double lambda, int n_x, int n_sp_x) {
  //calculate square root of P
  MatrixXd A = P_in.llt().matrixL();

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x, n_sp_x);

  //set first column of sigma point matrix
  Xsig.col(0)  = x_in;

  //set remaining sigma points
  for (int i = 0; i < n_x; i++)
  {
    Xsig.col(i+1)     = x_in + sqrt(lambda+n_x) * A.col(i);
    Xsig.col(i+1+n_x) = x_in - sqrt(lambda+n_x) * A.col(i);
  }

  //write result
  *Xsig_out = Xsig;
}

/**
 * @brief UKF::AugmentedSigmaPoints
 * @param Xsig_out Augmented sigma points
 * @param x state vector
 * @param P State covariance matrix
 * @param lambda
 * @param n_aug Lenght of augmented sigma point vector
 * @param n_sp_xaug Number of augmented sigmapoint vectors
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out, const VectorXd* x_in, const MatrixXd* P_in,
                               double lambda, int n_aug, int n_sp_xaug, double std_a, double std_yawdd) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug);
  x_aug.fill(0.0);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug, n_aug);
  P_aug.fill(0.0);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, n_sp_xaug);
  X_sig_aug.fill(0.0);

  //create augmented mean state
  x_aug.head(5) = x_in;

  //create augmented covariance matrix by inserting P_ matrix to top-left corner
  P_aug.topLeftCorner(5,5) = P_in;
  // Insert process noise to bottom right corner
  P_aug(5,5) = std_a*std_a;
  P_aug(6,6) = std_yawdd*std_yawdd;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;  // Mean state
  for (int i = 0; i< n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * L.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug;
}

/**
 * @brief UKF::SigmaPointPrediction
 * @param Xsig_out
 */
void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, MatrixXd* Xsig_aug, double delta_t, int n_sp_xaug) {
  float epsilon = 0.00000001;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, n_sp_xaug);

  // Variables to store predictions
  double px = 0;
  double py = 0;
  double v = 0;
  double psi = 0;
  double psi_dot = 0;

  //predict sigma points
  for(int i=0; i < n_sp_xaug; i++){
    // Variables to store (previous state) xk values for easier access to them
    float px_k = Xsig_aug(0, i);
    float py_k = Xsig_aug(1, i);
    float v_k = Xsig_aug(2, i);
    float psi_k = Xsig_aug(3, i);
    float psi_dot_k = Xsig_aug(4, i);
    float va_k = Xsig_aug(5, i);
    float psi_dotdot_k = Xsig_aug(6, i);

    //std::cout << "Avoid zero division..." << std::endl;
    //avoid division by zero (only for px, py and psi)
    if (Xsig_aug(4, i) > epsilon){
      px = px_k + (v_k / psi_dot_k) * ( sin(psi_k + (psi_dot_k * delta_t)) - sin(psi_k));  // predict px_k+1
      py = py_k + (v_k / psi_dot_k) * (-cos(psi_k + (psi_dot_k * delta_t)) + cos(psi_k));
      psi = psi_k + (psi_dot_k * delta_t);
    }
    else {
      // Predict px
      px = px_k + v_k * cos(psi_k) * delta_t;
      py = py_k + v_k * sin(psi_k) * delta_t;
      psi = psi_k;
    }
    // Predict remaining elements
    v = v_k;  // Predict v_k+1
    psi_dot = psi_dot_k;

    // Add noise
    px += 0.5 * pow(delta_t, 2) * cos(psi_k) * va_k;  // Add noise to px_k+1
    py += 0.5 * pow(delta_t, 2) * sin(psi_k) * va_k;
    v  += delta_t * va_k;
    psi += 0.5 * pow(delta_t, 2) * psi_dotdot_k;
    psi_dot += delta_t * psi_dotdot_k;

    //write predicted sigma points into right column
    Xsig_pred(0, i) = px;
    Xsig_pred(1, i) = py;
    Xsig_pred(2, i) = v;
    Xsig_pred(3, i) = psi;
    Xsig_pred(4, i) = psi_dot;
  }

  //write result
  *Xsig_out = Xsig_pred;

}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);

  //create vector for predicted state
  VectorXd x = VectorXd(n_x);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);


/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //set weights
  weights(0) = float(lambda) / float((lambda + n_aug));
  for (int i=1; i < (2*n_aug+1); i++){
    float w = 1. / (2 * (lambda + n_aug));
    weights(i) = w;
  }
  //predict state mean
  x.fill(0.0);
  for (int i=0; i < (2*n_aug+1); i++){
    x += Xsig_pred.col(i) * weights(i);
  }
  //predict state covariance matrix
  P.fill(0.0);
  for (int i=0; i < (2* n_aug + 1); i++){
      // Calculate difference
      VectorXd x_diff = Xsig_pred.col(i) - x;
      // Normalize angles
      while (x_diff(3) > M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;
      P = P + weights(i) * x_diff * x_diff.transpose();
  }

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
   double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  //radar measurement noise standard deviation radius in m
  double std_radr = 0.3;

  //radar measurement noise standard deviation angle in rad
  double std_radphi = 0.0175;

  //radar measurement noise standard deviation radius change in m/s
  double std_radrd = 0.1;

  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  for (int i=0; i < (n_z, 2 * n_aug + 1); i++){
      // Calculate p
      Zsig(0, i) = sqrt(pow(Xsig_pred(0,i), 2) + pow(Xsig_pred(1,i), 2));
      // Calculate phi and normalize angle
      Zsig(1, i) = atan(Xsig_pred(1,i) / Xsig_pred(0,i));
      //while (Zsig(1, i) > M_PI) Zsig(1, i)-=2.*M_PI;
      //while (Zsig(1, i) < M_PI) Zsig(1, i)+=2.*M_PI;
      // calculate p dot
      Zsig(2, i) = ((Xsig_pred(0,i) * cos(Xsig_pred(3,i)) * Xsig_pred(2,i)) + (Xsig_pred(1,i) * sin(Xsig_pred(3,i)) * Xsig_pred(2,i))) / Zsig(0, i);

  }
  //calculate mean predicted measurement
  //predict state mean
  z_pred.fill(0.0);
  for (int i=0; i < (2*n_aug+1); i++){
    z_pred += Zsig.col(i) * weights(i);
  }

  // R-matrix
  MatrixXd R = MatrixXd(3, 3);
  R << pow(std_radr, 2), 0, 0,
       0, pow(std_radphi, 2), 0,
       0, 0, pow(std_radrd, 2);

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i=0; i < (2* n_aug + 1); i++){
      // Calculate difference
      VectorXd z_diff = Zsig.col(i) - z_pred;
      // Normalize angles
      while (z_diff(1) > M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;
      S = S + weights(i) * z_diff * z_diff.transpose();
  }
  S += R;

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;

  //write result
  *z_out = z_pred;
  *S_out = S;
}

void UKF::UpdateState(VectorXd* x_out, MatrixXd* P_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
   double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  //create example matrix with predicted sigma points in state space
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create example vector for predicted state mean
  VectorXd x = VectorXd(n_x);
  x <<
     5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;

  //create example matrix for predicted state covariance
  MatrixXd P = MatrixXd(n_x,n_x);
  P <<
  0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
 -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
 -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

  //create example matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  Zsig <<
      6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
     0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
      2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;

  //create example vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred <<
      6.12155,
     0.245993,
      2.10313;

  //create example matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);
  S <<
      0.0946171, -0.000139448,   0.00407016,
   -0.000139448,  0.000617548, -0.000770652,
     0.00407016, -0.000770652,    0.0180917;

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<
      5.9214,   //rho in m
      0.2187,   //phi in rad
      2.0062;   //rho_dot in m/s

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  for (int i=0; i < (2 * n_aug +1); i++) {
    // Calculate difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();

  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // Residual
  VectorXd z_diff = z - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();


/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Updated state x: " << std::endl << x << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}

