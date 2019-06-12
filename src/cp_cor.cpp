/*--------Parallel correlation calculation----------*/
// Modified from https://systematicinvestor.github.io/Correlation-Rcpp
// [[Rcpp::depends(RcppParallel)]]
//
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

// pre-compute sum and stdev
struct cor_p1 :  public Worker {
    const RcppParallel::RMatrix<double> mat;
    const int rstart, rend, nperiod;

    RcppParallel::RVector<double> rsum, rstdev;

    cor_p1(const Rcpp::NumericMatrix& mat, const int rstart, const int rend,
            Rcpp::NumericVector rsum, Rcpp::NumericVector rstdev)
        : mat(mat), rstart(rstart), rend(rend), nperiod(rend - rstart),
        rsum(rsum), rstdev(rstdev) {  }

    void operator() (size_t begin, size_t end) {
        for (size_t c = begin; c < end; ++c) {
            double sum, sum2;
            sum = sum2 = 0;

            for (int r = rstart; r < rend; ++r) {
                double d = mat(r,c);
                sum += d;
                sum2 += pow(d,2);
            }
            rsum[c] = sum;
            rstdev[c] = sqrt(nperiod * sum2 - pow(sum,2));
        }
    }
};

// compute correlation
struct cor_p2 : public Worker {
    const RcppParallel::RMatrix<double> mat;
    const int rstart, rend, nperiod;
    const RcppParallel::RVector<double> sum, stdev;

    RcppParallel::RMatrix<double> rmat;

    cor_p2(const Rcpp::NumericMatrix& mat, const int rstart, const int rend,
            const Rcpp::NumericVector& sum, const Rcpp::NumericVector& stdev,
            Rcpp::NumericMatrix rmat)
        : mat(mat), rstart(rstart), rend(rend), nperiod(rend - rstart),
        sum(sum), stdev(stdev), rmat(rmat) {  }
    void operator()(size_t begin, size_t end) {
        for (size_t c1 = begin; c1 < end; c1++) {
            for (size_t c2 = 0; c2 < c1; c2++) {
                double sXY = 0;
                for (int r = rstart; r < rend; r++) {
                    sXY += mat(r, c1) * mat(r, c2);
                }

                double sX = sum[c1];
                double sY = sum[c2];
                double SDx = stdev[c1];
                double SDy = stdev[c2];
                rmat(c1, c2) = (nperiod * sXY - sX * sY) / (SDx * SDy);
                rmat(c2, c1) = rmat(c1, c2);
            }
        }
    }
};

Rcpp::NumericMatrix cp_cor_helper(const Rcpp::NumericMatrix& mat, const int rstart,
        const int rend) {
    int nc = mat.ncol();
    NumericVector rsum(nc), rstdev(nc);

    cor_p1 p1(mat, rstart, rend, rsum, rstdev);
    parallelFor(0, nc, p1);

    Rcpp::NumericMatrix rmat(nc, nc);

    cor_p2 p2(mat, rstart, rend, rsum, rstdev, rmat);
    parallelFor(0, nc, p2);

    return rmat;
}

// @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cp_cor(Rcpp::NumericMatrix mat) {
    return cp_cor_helper(mat, 0, mat.nrow());
}
