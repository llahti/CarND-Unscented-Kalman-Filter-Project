#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
}

VectorXd Tools::Polar2Cartesian(VectorXd &polar){
  VectorXd cart = VectorXd(4);
  cart << polar[0] * cos(polar[1]),
          polar[0] * sin(polar[1]),
          polar[2] * cos(polar[1]),
          polar[2] * sin(polar[1]);
  return cart;
}
