#include "WireCellUtil/LassoModel.h"
#include "WireCellUtil/ElasticNetModel.h"

#include "WireCellUtil/Eigen.h"
using namespace Eigen;

#include <iostream>
using namespace std;

void test_model(WireCell::LinearModel& m, MatrixXd& G, VectorXd& W);
void test_lasso(WireCell::LassoModel& m, MatrixXd& G, VectorXd& W);
void print_results(WireCell::LinearModel& m, VectorXd& C);

int main(int argc, char* argv[])
{
    // std::srand((unsigned int) time(0));

    const int N_CELL = 40;
    const int N_ZERO = 30;
    const int N_WIRE = int(N_CELL * 0.7);

    // initialize C vector: NCELL cells with NZERO zeros. (true charge in each cell)
    VectorXd C = VectorXd::Random(N_CELL) * 50 + VectorXd::Constant(N_CELL, 150);
    VectorXd r = N_CELL / 2 * (VectorXd::Random(N_ZERO) + VectorXd::Constant(N_ZERO, 1));
    for (int i = 0; i < N_ZERO; i++) {
        C(int(r(i))) = 0;
    }

    // initialize G matrix: N_WIRE rows and N_CELL columns. (geometry matrix)
    MatrixXd G = MatrixXd::Zero(N_WIRE, N_CELL);
    for (int i = 0; i < N_WIRE; i++) {
        VectorXd t = VectorXd::Random(N_CELL);
        for (int j = 0; j < N_CELL; j++) {
            G(i, j) = int(t(j) + 1);
        }
    }

    // W vector is the measured charge on wires.
    VectorXd W = G * C;

    cout << "geometry matrix:" << endl;
    cout << G << endl << endl;

    // cout << "measured charge on each wire:" << endl;
    // cout << W.transpose() << endl << endl;

    cout << "true charge of each cell:" << endl;
    cout << C.transpose() << endl << endl;
    double lambda = 0.01;

    // WireCell::ElasticNetModel m(lambda, 0.95, 100000, 1e-4);
    // test_model(m, G, W);
    // print_results(m, C);

    WireCell::LassoModel m2(lambda, 100000, 1e-3);
    test_lasso(m2, G, W);
    print_results(m2, C);

    return 0;
}

void test_model(WireCell::LinearModel& m, MatrixXd& G, VectorXd& W)
{
    m.SetData(G, W);
    m.Fit();
}

void test_lasso(WireCell::LassoModel& m, MatrixXd& G, VectorXd& W)
{
    m.SetData(G, W);

    // one can set the weight of each cell's regularization.
    // m.SetLambdaWeight(0, 10.);
    // m.SetLambdaWeight(2, 10.);
    // m.SetLambdaWeight(6, 10.);
    m.Fit();

    cout << "chi2: " << m.name << ": base=" << m.chi2_base() << ", l1=" << m.chi2_l1() << endl << endl;
}

void print_results(WireCell::LinearModel& m, VectorXd& C)
{
    VectorXd beta = m.Getbeta();

    cout << "fitted charge of each cell: " << m.name << endl;
    cout << beta.transpose() << endl << endl;

    // cout << "predicted charge on each wire: Lasso" << endl;
    // cout << m.Predict().transpose() << endl << endl;

    cout << "average residual charge difference per wire: " << m.name << ": " << m.MeanResidual() << endl;

    int nbeta = beta.size();

    int n_zero_true = 0;
    int n_zero_beta = 0;
    int n_zero_correct = 0;

    for (int i = 0; i < nbeta; i++) {
        if (fabs(C(i)) < 0.1) n_zero_true++;
        if (fabs(beta(i)) < 5) n_zero_beta++;
        if (fabs(C(i)) < 0.1 && (fabs(C(i) - beta(i)) < 10)) n_zero_correct++;
    }
    cout << "true zeros: " << n_zero_true << endl;
    cout << "fitted zeros: " << n_zero_beta << endl;
    cout << "correct fitted zeros: " << n_zero_correct << endl;

    cout << endl;
}
