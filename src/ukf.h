#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {

private:
  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;

public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* state covariance matrix for augmented sigmapoints
  MatrixXd P_aug_;

  ///* Measurement covariance matrix for radar measurements
  MatrixXd S_rad_;

  ///* sigma points matrix
  MatrixXd Xsig_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* augmented sigma points matrix
  MatrixXd Xsig_aug_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  /// Number of sigma-points for state matrix
  int n_sp_x_;

  /// Number of sigma-points for augmented state matrix
  int n_sp_xaug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Augmented Sigma point spreading parameter
  double lambda_aug_;

  ///* Length of radar measurement vector
  int n_z_rad_;

  ///* Length of lidar measurement vector
  int n_z_lidar_;

  ///* H matrix for Lidar update step
  MatrixXd H_lidar_;
  MatrixXd R_lidar_;  // Measurement Covariance matrix


  /**
   * Constructor
   */
  UKF();


  /**
   * Destructor
   */
  virtual ~UKF();

  void FirstUpdate(MeasurementPackage measurement_pack);


  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);


  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);


  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);


  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);


  /**
   * @brief GenerateSigmaPoints
   */
  void GenerateSigmaPoints();


  /**
   * @brief AugmentedSigmaPoints
   */
  void AugmentedSigmaPoints();

  void SigmaPointPrediction(double delta_t);

  void PredictMeanAndCovariance();

  void PredictLidarMeasurement(VectorXd* z_pred_out);

  void PredictRadarMeasurement(VectorXd* z_pred_out, MatrixXd* Zsig_out);

  void UpdateLidarState(MatrixXd* Zsig, VectorXd* z_pred, VectorXd* z);

  void UpdateRadarState(MatrixXd* Zsig, VectorXd* z_pred, VectorXd* z);
};

#endif /* UKF_H */
