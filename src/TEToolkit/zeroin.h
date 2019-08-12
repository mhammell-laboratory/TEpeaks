//
//  Prob.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/25/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____zeroin__
#define __TEToolkit_c____zeroin__

#include <stdio.h>
#include <vector>
#include <string>
#include "Constants.h"

extern "C" {
    double ppois_inv_lambda(double lambda,int k, double prob);
    double ppois(int k,double pval);
    
    double poisson_cdf ( unsigned int n, double lam, bool lower=false, bool log10=false );
    std::vector<double> convolve (std::vector<double> a, std::vector<double> b, std::string mode ="valid");
    
    int binomial_cdf_inv (double cdf,int a,double b);
    double standard_deviation(std::vector<int> vals, double avg);
    std::vector<double> correlate (std::vector<double> a, std::vector<double> b, int mode =1);
    void linspace(int start, int stop, int num,std::vector<int>& res);
    void smooth(std::vector<double>& x, int window_len=11, std::string window="hanning");

double zeroin2(			/* An estimate of the root */
               int kk,
               double prob,
               double ax,				/* Left border | of the range	*/
               double bx,				/* Right border| the root is seeked*/
               double fa, double fb,		/* f(a), f(b) */
               //double (*f)(double x, void *info),	/* Function under investigation	*/
               //void *info,				/* Add'l info passed on to f	*/
               double *Tol,			/* Acceptable tolerance		*/
               int *Maxit);

    double retrieve_lambda_bg(int k,double prob);
    
    //Signal
    std::vector<int> enforce_peakyness(std::vector<double> signal,std::vector<int> maximum);
    void internal_minima( std::vector<double> signal, std::vector<int> maximum,std::vector<int> &ret );
    bool is_valid_peak(std::vector<double> signal, int maximum);
    std::vector<double> hard_clip(std::vector<double> signal, int maximum);
    bool too_flat(std::vector<double> &signal);
    void maxima(std::vector<double> signal,std::vector<int> ret,int window_size=51);
    std::vector<double> savitzky_golay_order2(std::vector<double> signal,
                                              int window_size, int deriv=0);
    
    
}

#endif /* defined(__TEToolkit_c____zeroin__) */
