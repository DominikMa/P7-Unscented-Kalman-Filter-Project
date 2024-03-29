#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the inputs
  if (estimations.empty()) return rmse;
  if (estimations.size() != ground_truth.size()) return rmse;

  //accumulate squared residuals
  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd d = estimations[i] - ground_truth[i];
    rmse += d.cwiseProduct(d);
  }

  //calculate the mean
  rmse /= estimations.size();

  //calculate the squared root
  rmse = rmse.cwiseSqrt();

  //return the result
  return rmse;
}