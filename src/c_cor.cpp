// cor.cpp
//
// This script is part of MINA.
// Modified from http://systematicinvestor.github.io/Correlation-Rcpp
//
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// Define the type `asset_info`.
struct asset_info {
    double sum, sum2, stdev; // sX, sX2 and Sx'
};

// [corrlation matrix]
// (http://en.wikipedia.org/wiki/Correlation_and_dependence).
// n, sX, sY, sXY, sX2, sY2
// cor = ( n * sXY - sX * sY) / ( sqrt(n * sX2^2) * sqrt(n * sY2 - sY^2) )

inline asset_info compute_asset_info(const NumericMatrix& mat,
        const int icol, const int rstart, const int rend) {
    double sum, sum2;
    sum = sum2 = 0;

    // read in the column of the matrix and summing up
    for (int r = rstart; r < rend; r++) {
        double d = mat(r, icol);
        sum += d;
        sum2 += pow(d, 2);
    }

    // new object and assign the sum, sum2 and stdev values
    asset_info res;
        res.sum = sum;
        res.sum2 = sum2;
        res.stdev = sqrt((rend - rstart) * sum2 - pow(sum, 2));
    return res;
}

inline NumericMatrix c_cor_helper (const NumericMatrix& mat, const int rstart,
        const int rend) {

    int nc = mat.ncol();
    int nperiod = rend - rstart;

    // corretion between columns are calculated
    NumericMatrix rmat(nc, nc);

    vector<asset_info> info(nc);
    for (int c = 0; c < nc; c++) {
        info[c] = compute_asset_info(mat, c, rstart, rend);
    }
    for (int c1 = 0; c1 < nc; c1++) {
        for (int c2 = 0; c2 < c1; c2++) {
            double sXY = 0;
            for (int r = rstart; r < rend; r++){
                sXY += mat(r, c1) * mat(r, c2);
            }
            double sX = info[c1].sum;
            double sY = info[c2].sum;
            double SDx = info[c1].stdev;
            double SDy = info[c2].stdev;
            rmat(c1, c2) = (nperiod * sXY - sX * sY) / (SDx * SDy);
        }
    }
    return rmat;
 }

// [[Rcpp::export]]
NumericMatrix c_cor(NumericMatrix mat) {
    return c_cor_helper(mat, 0, mat.nrow());
}
