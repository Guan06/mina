/*-------------Parallel log-transformed variances calculation-----*/
// Modified from https://systematicinvestor.github.io/Correlation-Rcpp
//
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

// pre-compute log
struct log_trans : public Worker {
    const RcppParallel::RMatrix<double> mat;
    const int trans;
};
