#include <math.h>

//#define max(x, y) ((x)>(y)?(x):(y))
#define logspace_add(x, y) ((x)>(y)?(x):(y)) + log1p ( exp( -fabs(x - y) ) )
#define LOG10_E 0.43429448190325176

extern "C" {
double log10_poisson_cdf ( unsigned int n, double lam, short lower );
double log10_poisson_cdf_P_large_lambda ( unsigned int k, double lbd );
double log10_poisson_cdf_Q_large_lambda ( unsigned int k, double lbd );
double chi2_k1_cdf ( double x );
double chi2_k2_cdf ( double x );
double chi2_k4_cdf ( double x );
double log10_chi2_k1_cdf ( double x );
double log10_chi2_k2_cdf ( double x );
double log10_chi2_k4_cdf ( double x );
}
