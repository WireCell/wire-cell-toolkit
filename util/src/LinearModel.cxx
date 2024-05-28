#include "WireCellUtil/LinearModel.h"

#include "WireCellUtil/Eigen.h"
using namespace Eigen;

WireCell::LinearModel::LinearModel() {}

WireCell::LinearModel::~LinearModel() {}

VectorXd WireCell::LinearModel::Predict() { return _X * _beta; }

double WireCell::LinearModel::chi2_base() { return (_y - Predict()).squaredNorm(); }

double WireCell::LinearModel::MeanResidual() { return (_y - Predict()).norm() / _y.size(); }
