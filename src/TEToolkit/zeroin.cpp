/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999, 2001 the R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* from NETLIB c/brent.shar with max.iter, add'l info and convergence
   details hacked in by Peter Dalgaard */

/*************************************************************************
 *			    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,info,Tol,Maxit)
 *	double ax;			Root will be seeked for within
 *	double bx;			a range [ax,bx]
 *	double (*f)(double x, void *info); Name of the function whose zero
 *					will be seeked for
 *	void *info;			Add'l info passed to f
 *	double *Tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 *	int *Maxit;			Max. iterations
 *
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *	*Tol returns estimated precision
 *	*Maxit returns actual # of iterations, or -1 if maxit was
 *      reached without convergence.
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bisection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bisection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bisection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 *
 * NOTE:  uniroot() --> do_zeroin2()  --- in  ../main/optimize.c
 *					      ~~~~~~~~~~~~~~~~~~
 */

#include "zeroin.h"
#include <vector>
//#include <cfloat>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>




#include <boost/lexical_cast.hpp>

#include "myLog.h"
#include "cStatistics.h"


//using namespace std;

using boost::lexical_cast;


#define LSTEP  200
#define EXPTHRES exp(LSTEP)
#define EXPSTEP  exp(-LSTEP)

double sagv_coeff_51[51] ={ -2.26244344e-03,  -2.17194570e-03,  -2.08144796e-03,
    -1.99095023e-03,  -1.90045249e-03,  -1.80995475e-03,
    -1.71945701e-03,  -1.62895928e-03,  -1.53846154e-03,
    -1.44796380e-03,  -1.35746606e-03,  -1.26696833e-03,
    -1.17647059e-03,  -1.08597285e-03,  -9.95475113e-04,
    -9.04977376e-04,  -8.14479638e-04,  -7.23981900e-04,
    -6.33484163e-04,  -5.42986425e-04,  -4.52488688e-04,
    -3.61990950e-04,  -2.71493213e-04,  -1.80995475e-04,
    -9.04977376e-05,  -2.64877094e-19,   9.04977376e-05,
    1.80995475e-04,   2.71493213e-04,   3.61990950e-04,
    4.52488688e-04,   5.42986425e-04,   6.33484163e-04,
    7.23981900e-04,   8.14479638e-04,   9.04977376e-04,
    9.95475113e-04,   1.08597285e-03,   1.17647059e-03,
    1.26696833e-03,   1.35746606e-03,   1.44796380e-03,
    1.53846154e-03,   1.62895928e-03,   1.71945701e-03,
    1.80995475e-03,   1.90045249e-03,   1.99095023e-03,
    2.08144796e-03,   2.17194570e-03,   2.26244344e-03};

double sagv_coeff_201[201] = { -1.47775972e-04,  -1.46298212e-04,  -1.44820452e-04,
    -1.43342692e-04,  -1.41864933e-04,  -1.40387173e-04,
    -1.38909413e-04,  -1.37431654e-04,  -1.35953894e-04,
    -1.34476134e-04,  -1.32998374e-04,  -1.31520615e-04,
    -1.30042855e-04,  -1.28565095e-04,  -1.27087336e-04,
    -1.25609576e-04,  -1.24131816e-04,  -1.22654056e-04,
    -1.21176297e-04,  -1.19698537e-04,  -1.18220777e-04,
    -1.16743018e-04,  -1.15265258e-04,  -1.13787498e-04,
    -1.12309738e-04,  -1.10831979e-04,  -1.09354219e-04,
    -1.07876459e-04,  -1.06398700e-04,  -1.04920940e-04,
    -1.03443180e-04,  -1.01965420e-04,  -1.00487661e-04,
    -9.90099010e-05,  -9.75321413e-05,  -9.60543816e-05,
    -9.45766218e-05,  -9.30988621e-05,  -9.16211024e-05,
    -9.01433427e-05,  -8.86655830e-05,  -8.71878233e-05,
    -8.57100635e-05,  -8.42323038e-05,  -8.27545441e-05,
    -8.12767844e-05,  -7.97990247e-05,  -7.83212650e-05,
    -7.68435052e-05,  -7.53657455e-05,  -7.38879858e-05,
    -7.24102261e-05,  -7.09324664e-05,  -6.94547067e-05,
    -6.79769469e-05,  -6.64991872e-05,  -6.50214275e-05,
    -6.35436678e-05,  -6.20659081e-05,  -6.05881484e-05,
    -5.91103887e-05,  -5.76326289e-05,  -5.61548692e-05,
    -5.46771095e-05,  -5.31993498e-05,  -5.17215901e-05,
    -5.02438304e-05,  -4.87660706e-05,  -4.72883109e-05,
    -4.58105512e-05,  -4.43327915e-05,  -4.28550318e-05,
    -4.13772721e-05,  -3.98995123e-05,  -3.84217526e-05,
    -3.69439929e-05,  -3.54662332e-05,  -3.39884735e-05,
    -3.25107138e-05,  -3.10329540e-05,  -2.95551943e-05,
    -2.80774346e-05,  -2.65996749e-05,  -2.51219152e-05,
    -2.36441555e-05,  -2.21663957e-05,  -2.06886360e-05,
    -1.92108763e-05,  -1.77331166e-05,  -1.62553569e-05,
    -1.47775972e-05,  -1.32998374e-05,  -1.18220777e-05,
    -1.03443180e-05,  -8.86655830e-06,  -7.38879858e-06,
    -5.91103887e-06,  -4.43327915e-06,  -2.95551943e-06,
    -1.47775972e-06,   2.11874408e-20,   1.47775972e-06,
    2.95551943e-06,   4.43327915e-06,   5.91103887e-06,
    7.38879858e-06,   8.86655830e-06,   1.03443180e-05,
    1.18220777e-05,   1.32998374e-05,   1.47775972e-05,
    1.62553569e-05,   1.77331166e-05,   1.92108763e-05,
    2.06886360e-05,   2.21663957e-05,   2.36441555e-05,
    2.51219152e-05,   2.65996749e-05,   2.80774346e-05,
    2.95551943e-05,   3.10329540e-05,   3.25107138e-05,
    3.39884735e-05,   3.54662332e-05,   3.69439929e-05,
    3.84217526e-05,   3.98995123e-05,   4.13772721e-05,
    4.28550318e-05,   4.43327915e-05,   4.58105512e-05,
    4.72883109e-05,   4.87660706e-05,   5.02438304e-05,
    5.17215901e-05,   5.31993498e-05,   5.46771095e-05,
    5.61548692e-05,   5.76326289e-05,   5.91103887e-05,
    6.05881484e-05,   6.20659081e-05,   6.35436678e-05,
    6.50214275e-05,   6.64991872e-05,   6.79769469e-05,
    6.94547067e-05,   7.09324664e-05,   7.24102261e-05,
    7.38879858e-05,   7.53657455e-05,   7.68435052e-05,
    7.83212650e-05,   7.97990247e-05,   8.12767844e-05,
    8.27545441e-05,   8.42323038e-05,   8.57100635e-05,
    8.71878233e-05,   8.86655830e-05,   9.01433427e-05,
    9.16211024e-05,   9.30988621e-05,   9.45766218e-05,
    9.60543816e-05,   9.75321413e-05,   9.90099010e-05,
    1.00487661e-04,   1.01965420e-04,   1.03443180e-04,
    1.04920940e-04,   1.06398700e-04,   1.07876459e-04,
    1.09354219e-04,   1.10831979e-04,   1.12309738e-04,
    1.13787498e-04,   1.15265258e-04,   1.16743018e-04,
    1.18220777e-04,   1.19698537e-04,   1.21176297e-04,
    1.22654056e-04,   1.24131816e-04,   1.25609576e-04,
    1.27087336e-04,   1.28565095e-04,   1.30042855e-04,
    1.31520615e-04,   1.32998374e-04,   1.34476134e-04,
    1.35953894e-04,   1.37431654e-04,   1.38909413e-04,
    1.40387173e-04,   1.41864933e-04,   1.43342692e-04,
    1.44820452e-04,   1.46298212e-04,   1.47775972e-04};

double sagv_coeff_101[101] = {
    -7.38879858e-05,
    -7.24102261e-05,  -7.09324664e-05,  -6.94547067e-05,
    -6.79769469e-05,  -6.64991872e-05,  -6.50214275e-05,
    -6.35436678e-05,  -6.20659081e-05,  -6.05881484e-05,
    -5.91103887e-05,  -5.76326289e-05,  -5.61548692e-05,
    -5.46771095e-05,  -5.31993498e-05,  -5.17215901e-05,
    -5.02438304e-05,  -4.87660706e-05,  -4.72883109e-05,
    -4.58105512e-05,  -4.43327915e-05,  -4.28550318e-05,
    -4.13772721e-05,  -3.98995123e-05,  -3.84217526e-05,
    -3.69439929e-05,  -3.54662332e-05,  -3.39884735e-05,
    -3.25107138e-05,  -3.10329540e-05,  -2.95551943e-05,
    -2.80774346e-05,  -2.65996749e-05,  -2.51219152e-05,
    -2.36441555e-05,  -2.21663957e-05,  -2.06886360e-05,
    -1.92108763e-05,  -1.77331166e-05,  -1.62553569e-05,
    -1.47775972e-05,  -1.32998374e-05,  -1.18220777e-05,
    -1.03443180e-05,  -8.86655830e-06,  -7.38879858e-06,
    -5.91103887e-06,  -4.43327915e-06,  -2.95551943e-06,
    -1.47775972e-06,   2.11874408e-20,   1.47775972e-06,
    2.95551943e-06,   4.43327915e-06,   5.91103887e-06,
    7.38879858e-06,   8.86655830e-06,   1.03443180e-05,
    1.18220777e-05,   1.32998374e-05,   1.47775972e-05,
    1.62553569e-05,   1.77331166e-05,   1.92108763e-05,
    2.06886360e-05,   2.21663957e-05,   2.36441555e-05,
    2.51219152e-05,   2.65996749e-05,   2.80774346e-05,
    2.95551943e-05,   3.10329540e-05,   3.25107138e-05,
    3.39884735e-05,   3.54662332e-05,   3.69439929e-05,
    3.84217526e-05,   3.98995123e-05,   4.13772721e-05,
    4.28550318e-05,   4.43327915e-05,   4.58105512e-05,
    4.72883109e-05,   4.87660706e-05,   5.02438304e-05,
    5.17215901e-05,   5.31993498e-05,   5.46771095e-05,
    5.61548692e-05,   5.76326289e-05,   5.91103887e-05,
    6.05881484e-05,   6.20659081e-05,   6.35436678e-05,
    6.50214275e-05,   6.64991872e-05,   6.79769469e-05,
    6.94547067e-05,   7.09324664e-05,   7.24102261e-05,
    7.38879858e-05};


double binomial_pdf( int x, int a, double b )
{
    //"""binomial PDF by H. Gene Shin
    
    if (a<1){
        return 0;
    }
    if( x < 0 or a < x)
    {
        return 0;
    }
    if( b==0)
    {
        if (x==0)
        {
            return 1;
        }
        else{
            return 0;
        }
    }
    if (b==1){
        if (x==a)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }

    double p;
    int mn;
    long mx;
//    int q;

    if( x > a - x)
    {
        p = 1.0 - b;
        mn = a-x;
        mx = x;
    }
    else
    {
        p=b;
        mn=x;
        mx=a-x;
    }
    double pdf = 1.0;
    int t = 0;
    
    for(int q =1; q < mn+1; q++)
    {
        pdf *= (a-q+1.0) * p/(mn-q+1.0);
        if (pdf < 1e-100)
        {
            while (pdf < 1e-3)
            {
                pdf /= 1.0-p;
                t -= 1;
            }
        }
        if (pdf > 1e+100)
        {
            while (pdf > 1e+3 and t<mx)
            {
                pdf *= 1.0-p;
                t+=1;
            }
        }
    }
    for(int  i=0; i < mx -t; i++ )
    {
        pdf *= 1.0-p;
    }
    
    //pdf=float("%.10e" % pdf);
    
    return pdf;
}

int binomial_cdf_inv ( double cdf, int a, double b )
{
    //"""BINOMIAL_CDF_INV inverts the binomial CDF. For lower tail only!

    
    if (cdf < 0 or cdf >1)
    {
        error("CDF must >= 0 or <= 1\n");
        std::exit(1);
    }

    double cdf2 = 0.0;
    

    for(int i = 0; i < a+1; i++)
    {
        double pdf = binomial_pdf (i,a,b);
        cdf2 = cdf2 + pdf;
        if (cdf < cdf2)
        {
            return i;
        }
    }
    return a;
}

// Modified from http://www.scipy.org/Cookbook/SavitzkyGolay
// positive window_size not enforced anymore
// needs sane input paramters, window size > 4
// switched to double precision for internal accuracy
std::vector<double> savitzky_golay_order2(std::vector<double> signal,
                     int window_size, int deriv)
{
    /*"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    */
    
    //np.ndarray[np.int64_t, ndim=2] b;
    //pad the signal at the extremes with values taken from the signal itself
    std::vector<double> firstvals, lastvals, ret;
    
    //np.ndarray[np.float64_t, ndim=1] m;
    
    if (window_size % 2 != 1){
        window_size += 1;
    }
    
    int half_window = (int)(window_size - 1) / 2;
    
    std::vector<double> m;
    
    if (window_size == 51) {
        
        for (int k = 50; k >= 0; k++) {
            m.push_back(sagv_coeff_51[k]);
        }
    }
    if (window_size == 101) {
        for (int k = 100; k >= 0; k--) {
            m.push_back(sagv_coeff_101[k]);
        }
    }
    
    //m = np.linalg.pinv(b).A[deriv];
    
    // pad the signal at the extremes with values taken from the signal itself
    //std::vector<double> firstvals;
    //std::vector<double> lastvals;
    
    for (size_t k=half_window; k > 0; k --) {
        firstvals.push_back(signal[0] - std::abs(signal[k] - signal[0]));
    }
    for (int k = -1; k > -1 * half_window -1 ; k--) {
        lastvals.push_back(signal[signal.size()-1] + std::abs(signal[signal.size()+k] - signal[signal.size()-1]));
    }
    //firstvals = signal[0] - np.abs(signal[1:half_window+1][::-1] - signal[0]);
    //lastvals = signal[-1] + np.abs(signal[-half_window-1:-1][::-1] - signal[-1]);
    
    //signal = np.concatenate((firstvals, signal, lastvals));
    //ret = np.convolve( m[::-1], signal.astype('float64'), mode='valid').astype('float32');
    
    std::vector<double> ss;
    
    ss.insert(ss.end(),firstvals.begin(),firstvals.end());
    ss.insert(ss.end(),signal.begin(),signal.end());
    ss.insert(ss.end(),lastvals.begin(),lastvals.end());
    
    ret = convolve(m, ss, "valid");
    
    return ret ;
}


void maxima(std::vector<double> signal,std::vector<int> ret,int window_size)
{
    /*"""return the local maxima in a signal after applying a 2nd order
    Savitsky-Golay (polynomial) filter using window_size specified  
    """*/

    window_size = window_size/2*2+1; // to make a even number

    std::vector<double> smoothed = savitzky_golay_order2(signal, window_size, 1);
    
    for (size_t j=0; j < signal.size()-1; j++) {
        if (smoothed[j+1] < 0 && smoothed[j] >=0) {
            ret.push_back(j);
        }
    }
    //ret = np.where(np.diff(np.sign()) <= -1)[0].astype('int32');
}

// require at least 6 different float values -- prevents broad flat peaks
bool too_flat(std::vector<double> &signal)
{
//    """return whether signal has at least 6 unique values
    std::vector<double>::iterator it;
    it = std::unique (signal.begin(), signal.end());
    
    signal.resize( std::distance(signal.begin(),it) );
    
    if (signal.size() > 6) {
        return false;
    }
    else {
        return true;
    }
    
}

// hard clip a region with negative values
std::vector<double> hard_clip(std::vector<double> signal, int maximum)
{
//    """clip the signal in both directions at the nearest values <= 0 to position maximum
    
    int left = 0;
    int right = signal.size();
    
    // clip left
    for(int i=maximum -1; i >=0; i --)
    {
        if (signal[i] < 0)
        {
            left = i;
            break;
        }
    }
    for(int i = maximum; i < right; i++ )
    {
        if (signal[i] < 0)
        {
            right = i;
            break ;
        }
    }
    std::vector<double> s(signal.begin() + left, signal.begin() + right) ;
    
    return s;
}

// hardcoded minimum peak width = 50
bool is_valid_peak(std::vector<double> signal, int maximum)
{
    std::vector<double> s = hard_clip(signal, maximum);
    
    if (s.size() < 50)
    {
        return false;
    }
    else {
        if (too_flat(s)) {
            return false;
        }
    }
    return true;

}

void internal_minima( std::vector<double> signal, std::vector<int> maximum,std::vector<int> &ret )
{

    int n = maximum.size();
    //int i, v, v2;
    
    if (n > 1)
    {
        
        //ret = np.zeros(n - 1, 'int32');
        int pos1 = maximum[0];
        for (int i=0; i < n-1; i++)
        {
            int pos2 = maximum[i + 1];
            double m = *std::min_element(signal.begin()+pos1,signal.begin()+pos2);
            for (int j = pos1; j < pos2; j++) {
                if (signal[j] == m) {
                    ret.push_back(j);
                    break;
                }
            }
            //ret.push_back(np.where(signal[pos1:pos2] == signal[pos1:pos2].min())[0][0] + pos1);
            pos1 = pos2;
        }
        
    }
    
}

double __poisson_cdf (int k,double a)
{
    //Poisson CDF For small lambda. If a > 745, this will return
    //incorrect result.
    
    if (k < 0){
        return 0 ;                       // special cases
    }
    double next = exp( -a );
    double cdf = next;
    
    for(int i=1;i < k+1; i++){
        double last = next;
        next = last * a / i;
        cdf = cdf + next ;
    }
    if (cdf > 1){
        return 1;
    }
    else {
        return cdf;
    }
}

double __poisson_cdf_large_lambda (int k, double a )
{
    //Slower poisson cdf for large lambda. ( > 700 )
    assert(a > 700);
    
    if (k < 0){
        return 0.0;
    }
    
    if (k < 0){
        return 0 ;                      // special cases
    }
    int num_parts = int(a/LSTEP);
    
    double last_part = int(a) % LSTEP ;
    
    double lastexp = std::exp(-last_part) ;
    double next = EXPSTEP ;
    num_parts -= 1 ;
    double cdf = next;
    
    for(int i=1;i<k+1; i++)
    {
        double last = next ;
        next = last * a / i ;
        cdf = cdf + next;
        
        if (next > EXPTHRES or cdf > EXPTHRES)
        {
            if(num_parts>=1){
                cdf *= EXPSTEP;
                next *= EXPSTEP;
                num_parts -= 1;
            }
            else{
                cdf *= lastexp ;
                lastexp = 1;
            }
        }
    }
    
    for(int i=0 ;i < num_parts;i++){
        cdf *= EXPSTEP ;
    }
    cdf *= lastexp;
    return cdf;
}

double __poisson_cdf_Q ( unsigned int k, double a )
{
    /*"""internal Poisson CDF evaluater for upper tail with small
     lambda.
     
     Parameters:
     k	: observation
     a	: lambda
     """*/
    unsigned int i;
    
    if (k < 0){
        return 1.0 ;                       // special cases
    }
    double nextcdf ;
    
    nextcdf = exp( -1 * a );
    double lastcdf;
    
    for (i=1; i <= k ; i++)
    {
        lastcdf = nextcdf;
        nextcdf = lastcdf * a / i;
    }
    
    double cdf = 0.0;
    i = k+1;
    while (nextcdf >0.0){
        lastcdf = nextcdf;
        nextcdf = lastcdf * a / i;
        cdf += nextcdf;
        i+=1;
    }
    return cdf;
}

double __poisson_cdf_Q_large_lambda ( unsigned int k, double a )
{
    /* """Slower internal Poisson CDF evaluater for upper tail with large
     lambda.
     
     Parameters:
     k	: observation
     a	: lambda
     """*/
    assert(a > 700);
    
    if (a <= 700) {
        error("This function works for lambda > 700.");
        std::exit(1);
    }
    
    if (k < 0){
        return 1.0 ;                       // special cases
    }
    unsigned int num_parts = int(a/LSTEP);
    double lastexp = exp(-1 * (((int)a) % LSTEP) );
    
    double nextcdf = EXPSTEP;
    unsigned int i;
    double lastcdf;
    
    num_parts -= 1;
    
    for( i=1; i <= k ; i++)
    {
        lastcdf = nextcdf;
        nextcdf = lastcdf * a / i;
        if (nextcdf > EXPTHRES){
            if (num_parts>=1)
            {
                nextcdf *= EXPSTEP;
                num_parts -= 1;
            }
            else{
                // simply raise an error
                error("Unexpected error in __poisson_cdf_Q_large_lambda.");
                std::exit(1);
                //cdf *= lastexp
                //lastexp = 1
            }
        }
    }
    
    double cdf = 0.0;
    i = k+1;
    while (nextcdf >0.0){
        lastcdf = nextcdf;
        nextcdf = lastcdf * a / i;
        cdf += nextcdf;
        i+=1;
        
        if (nextcdf > EXPTHRES || cdf > EXPTHRES){
            if (num_parts>=1){
                cdf *= EXPSTEP;
                nextcdf *= EXPSTEP;
                num_parts -= 1;
            }
            else{
                cdf *= lastexp;
                lastexp = 1;
            }
        }
    }
    
    for(i=0;i < num_parts;i++)
    {
        cdf *= EXPSTEP;
    }
    cdf *= lastexp;
    
    return cdf;
}




std::vector<int> enforce_peakyness(std::vector<double> signal,std::vector<int> maximum)
{
    /*"""requires peaks described by a signal and a set of points where the signal
    is at a maximum to meet a certain set of criteria
    
    maximum which do not meet the required criteria are discarded
    
    criteria:
        for each peak:
            calculate a threshold of the maximum of its adjacent two minima
                plus the sqrt of that value
            subtract the threshold from the region bounded by those minima
            clip that region if negative values occur inside it
            require it be > 50 bp in width
            require that it not be too flat (< 6 unique values) 
    """*/
    
    std::vector<int> minima;
    
    internal_minima(signal, maximum, minima);
    
    std::vector<double> new_signal;
    
    int n = minima.size();
    //float threshold;
    std::vector<int> peaky_maxima(maximum.begin(),maximum.end());
    int j = 0;
    
    if (n == 0){
        return maximum;
    }

    double threshold = signal[minima[0]];
    
    threshold += sqrt(threshold);
    
    for (int i = 0; i < minima[0]; i++ )
    {
        new_signal.push_back(signal[i] - threshold - sqrt(threshold));
    }
    
    if (is_valid_peak(new_signal, maximum[0]))
    {
        peaky_maxima[0] = maximum[0];
        j += 1;
    }
    
    for (int i=0;i < n-1; i ++)
    {
        threshold = std::max(signal[minima[i]], signal[minima[i + 1]]);
        threshold += sqrt(threshold);
        new_signal.clear();
        for (int k = minima[i]; k < minima[i+1]; k++) {
            new_signal.push_back(signal[k] - threshold);
        }
        
        int new_maximum = maximum[i+1] - minima[i];
        if (is_valid_peak(new_signal, new_maximum))
        {
            peaky_maxima[j] = maximum[i + 1];
            j += 1;
        }
    }
    threshold =  signal[minima[-1]];
    threshold += sqrt(threshold);
    for (size_t k = minima[minima.size()-1]; k < signal.size(); k++) {
        new_signal.push_back(signal[k] - threshold);
    }
    
    int new_maximum = maximum[-1] - minima[-1];
    
    if (is_valid_peak(new_signal, new_maximum))
    {
        peaky_maxima[j] = maximum[-1];
        j += 1;
    }
    peaky_maxima.resize(j);
    
    return peaky_maxima;
}



double standard_deviation(std::vector<int> vals, double avg)
{
    double sum =0 ;
    
    for (size_t i=0; i < vals.size(); i++) {
        
        sum += std::pow((vals[i] - avg),2);
    }
    return std::sqrt(sum)/(vals.size()-1);
}

void linspace(int start, int stop, int num,std::vector<int>& res)
{
    
    int step = (int)(1.0 * (stop - start +1)/num);
    
    res.push_back(start);
    for (int i= 1; i < num -1 ; i++) {
        res.push_back(start + step *i);
    }
    res.push_back(stop);
    
    
}

void shift2(std::vector<double> d, int s,std::vector<double>& c) {
    //std::vector<double> c(d.size(),0);
    
    for (size_t i = 0; i < c.size(); i++) {
        if ( s - i >= 0 )
            c[i] = d[s - i];
    }
    
}

void shift(std::vector<double> d, int s,std::vector<double>& c) {
    //std::vector<double> c(d.size(),0);
    
    for (size_t i = 0; i < d.size(); i++) {
        if ((s + i >= 0)&&(s+i < d.size()))
            c[i] = d[s + i];
    }
    
}

int dot(std::vector<double> a,std::vector<double> b)
{
    size_t aLength = a.size();
    
    if (aLength != b.size()) {
        std::cout << "ERROR: Vectors must be of equal length in dot product." << std::endl;
        return 0;
    }
    double sum = 0;
    for (size_t i = 0; i < aLength; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}


std::vector<double> convolve(std::vector<double> a, std::vector<double> b, std::string mode )
{
    /* for(auto v : a){
     std::cout << v << "\t";
     }
     std::cout << "\n";
     for(auto v : b){
     std::cout << v << "\t";
     }
     std::cout << "\n";*/
    
    //int i = -1 * b.size() + 1;
    
    //std::cout << i << std::endl;
    std::vector<double> res;
    
    int max_size,min_size;
    
    if (a.size() >= b.size()) {
        max_size = a.size();
        min_size = b.size();
        
        for ( int i = min_size - 1 ; i < max_size ; i++ ){
            std::cout << "i = " << i << std::endl;
            
            std::vector<double> slidingY(min_size,0);
            shift2(a, i,slidingY);
            
            /*  for(auto v : slidingY){
             std::cout << v << "\t";
             }
             std::cout << "\n";
             std::cout << "dot: " << dot(b,slidingY) << std::endl;*/
            
            res.push_back( dot(b,slidingY));
            //i += 1;
        }
    }
    if (a.size() < b.size() ) {
        max_size = b.size();
        min_size = a.size();
        
        for ( int i = min_size -1; i <  max_size ; i++ ){
            std::cout << "i = " << i << std::endl;
            
            std::vector<double> slidingY(min_size,0);
            shift2(b, i,slidingY);
            
            /* for(auto v : slidingY){
             std::cout << v << "\t";
             }
             std::cout << "\n";
             std::cout << "dot: " << dot(a,slidingY) << std::endl;
             //i += 1;*/
            res.push_back( dot(b,slidingY));
        }
    }
    return res;
}



std::vector<double> correlate (std::vector<double> a, std::vector<double> b, int mode )
{
    /*for(auto v : a){
     std::cout << v << "\t";
     }
     std::cout << "\n";
     for(auto v : b){
     std::cout << v << "\t";
     }
     std::cout << "\n";*/
    
    std::vector<double> cor;
    
    int size = b.size();
    
    for ( int i = size -1; i >  -1 * size; i-- ){
        //std::cout << "i = " << i << std::endl;
        
        std::vector<double> slidingY(b.size(),0);
        shift(b, i,slidingY);
        
        /* for(auto v : slidingY){
         std::cout << v << "\t";
         }
         std::cout << "\n";*/
        //std::cout << "dot: " << dot(a,slidingY) << std::endl;
        
        cor.push_back(dot(a,slidingY));
        
    }
    return cor;
}

// smooth function from SciPy cookbook: http://www.scipy.org/Cookbook/SignalSmooth

void smooth(std::vector<double>& x, int window_len, std::string window)
{
    /* smooth the data using a window with requested size.
     
     This method is based on the convolution of a scaled window with the signal.
     The signal is prepared by introducing reflected copies of the signal
     (with the window size) in both ends so that transient parts are minimized
     in the begining and end part of the output signal.
     input:
     x: the input signal
     window_len: the dimension of the smoothing window; should be an odd integer
     window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
     flat window will produce a moving average smoothing.
     
     output:
     the smoothed signal
     
     example:
     
     t=linspace(-2,2,0.1)
     x=sin(t)+randn(len(t))*0.1
     y=smooth(x)
     
     see also:
     
     numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
     scipy.signal.lfilter
     
     TODO: the window parameter could be the window itself if an array instead of a string
     NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
     */
    
    if( (int)x.size() < window_len){
        std::cout << "Input vector needs to be bigger than window size." << std::endl;
        std::exit(1);
    }
    if (window_len < 3){
        return ;
    }
    if( window != "flat" || window != "hanning" || window != "hamming" || window != "bartlett" || window != "blackman"){
        std::cout << "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'" << std::endl;
        std::exit(1);
    }
    std::vector<double> s;
    for (int i= window_len-1; i > 0; i--) {
        s.push_back(x[i]);
    }
    for (size_t i=0; i < x.size(); i ++) {
        s.push_back(x[i]);
    }
    for (size_t i = x.size()-1; i > x.size() - window_len ; i--) {
        s.push_back(x[i]);
    }
    //=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    
    std::vector<double> w;
    double w_sum = 0.0;
    if (window == "flat") //moving average
    {   for(int i = 0; i < window_len ; i ++)
    {
        w.push_back(1);
        w_sum += 1.0;
    }
    }
    else{
        //w=eval('np.'+window+'(window_len)');
        std::cout << "Window is not 'flat'" << std::endl;
        std::exit(1);
    }
    
    for(int i = 0; i < window_len ; i ++)
    {
        w[i] = w[i] /w_sum;
    }
    
    std::vector<double> y = convolve(s,w);
    
    x.clear();
    
    for (int i= window_len/2; i < window_len - window_len/2 ; i++)
    {
        x.push_back(y[i]);
    }
}


double poisson_cdf ( unsigned int n, double lam, bool lower, bool log10 )
{
    /*"""Poisson CDF evaluater.

    This is a more stable CDF function. It can tolerate large lambda
    value. While the lambda is larger than 700, the function will be a
    little slower.

    Parameters:
    n     : your observation
    lam   : lambda of poisson distribution
    lower : if lower is False, calculate the upper tail CDF, otherwise, to calculate lower tail; Default is False.
    log10 : if log10 is True, calculation will be in log space. Default is False.
    """*/
    
    //assert lam > 0.0, "Lambda must > 0, however we got %d" % lam
    if (lam <= 0) {
        error("Lambda must > 0, however we got " + lexical_cast<std::string> (lam));
        std::exit(1);
    }

    if (log10){
        if (lower){
            // lower tail
            return log10_poisson_cdf_P_large_lambda(n, lam);
        }
        else{
            // upper tail
            return log10_poisson_cdf_Q_large_lambda(n, lam);
        }
    }
    
    if (lower){
        if (lam > 700){
            return __poisson_cdf_large_lambda (n, lam);
        }
        else{
            return __poisson_cdf(n,lam);
        }
    }
    else{
        // upper tail
        if (lam > 700){
            return __poisson_cdf_Q_large_lambda (n, lam);
        }
        else {
            return __poisson_cdf_Q(n,lam);
        }
    }
}


double ppois(int x, double lam)
{
    /*Poisson CDF evaluater. This is a more stable CDF function. It can tolerate large lambda
     value. While the lambda is larger than 700, the function will be a
     little slower.
     
     Parameters:
     n     : your observation
     lam   : lambda of poisson distribution
     lower : if lower is False, calculate the upper tail CDF
     */
    int k = int(x);
    
    //if (lam <= 0.0){
    //    std::cout << "Lambda must > 0" << std::endl;
    //}
    
    bool lower = true;
    if (lower){
        if (lam > 700){
            return __poisson_cdf_large_lambda (k, lam);
        }
        else{
            return __poisson_cdf(k,lam);
        }
    }
    else {
        return -1;
    }
    //else{
    //   if lam > 700:
    //      return __poisson_cdf_Q_large_lambda (k, lam)
    //  else:
    //     return __poisson_cdf_Q(k,lam)
    //}
    
}


double ppois_inv_lambda(double lambda,int k, double prob)
{
    return ppois(k, lambda) - prob;
}

/* zeroin2() is faster for "expensive" f(), in those typical cases where
 *             f(ax) and f(bx) are available anyway : */

double zeroin2(			/* An estimate of the root */
               int kk, /* number of occurrence */
               double prob, /* probability */
               double ax,				/* Left border | of the range	*/
               double bx,				/* Right border| the root is seeked*/
               double fa, double fb,		/* f(a), f(b) */
               // double (*f)(double x, void *info),	/* Function under investigation	*/
               //void *info,				/* Add'l info passed on to f	*/
               double *Tol,			/* Acceptable tolerance	=2.220446e-16	*/
               int *Maxit)				/* Max # of iterations =100*/
{
    double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
    double tol;
    int maxit;
    
    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = * Tol;
    
    /* First test if we have found a root at an endpoint */
    if(fa == 0.0) {
        *Tol = 0.0;
        *Maxit = 0;
        return a;
    }
    if(fb ==  0.0) {
        *Tol = 0.0;
        *Maxit = 0;
        return b;
    }
    
    while(maxit--)		/* Main iteration loop	*/
    {
        double prev_step = b-a;		/* Distance from the last but one
                                     to the last approximation	*/
        double tol_act;			/* Actual tolerance		*/
        double p;			/* Interpolation step is calcu- */
        double q;			/* lated in the form p/q; divi-
                             * sion operations is delayed
                             * until the last moment	*/
        double new_step;		/* Step at this iteration	*/
        
        if( std::fabs(fc) < std::fabs(fb) )
        {				/* Swap data for b to be the	*/
            a = b;  b = c;  c = a;	/* best approximation		*/
            fa=fb;  fb=fc;  fc=fa;
        }
        tol_act = 2*EPSILON*std::fabs(b) + tol/2;
        new_step = (c-b)/2;
        
        if( std::fabs(new_step) <= tol_act || fb == (double)0 )
        {
            *Maxit -= maxit;
            *Tol = std::fabs(c-b);
            return b;			/* Acceptable approx. is found	*/
        }
        
        /* Decide if the interpolation can be tried	*/
        if( std::fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
           && std::fabs(fa) > std::fabs(fb) ) {	/* and was in true direction,
                                                 * Interpolation may be tried	*/
            register double t1,cb,t2;
            cb = c-b;
            if( a==c ) {		/* If we have only two distinct	*/
                /* points linear interpolation	*/
                t1 = fb/fa;		/* can only be applied		*/
                p = cb*t1;
                q = 1.0 - t1;
            }
            else {			/* Quadric inverse interpolation*/
                
                q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
                p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
                q = (q-1.0) * (t1-1.0) * (t2-1.0);
            }
            if( p >(double)0 )		/* p was calculated with the */
                q = -q;			/* opposite sign; make p positive */
            else			/* and assign possible minus to	*/
                p = -p;			/* q				*/
            
            if( p < (0.75*cb*q-std::fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
               && p < std::fabs(prev_step*q/2) )	/* and isn't too large	*/
                new_step = p/q;			/* it is accepted
                                         * If p/q is too large then the
                                         * bisection procedure can
                                         * reduce [b,c] range to more
                                         * extent */
        }
        
        if( std::fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
            if( new_step > (double)0 )	/* than tolerance		*/
                new_step = tol_act;
            else
                new_step = -tol_act;
        }
        a = b;	fa = fb;			/* Save the previous approx. */
        b += new_step;
        fb = ppois_inv_lambda(b,kk,prob);	/* Do step to a new approxim. */
        
        if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
            /* Adjust c for it to have a sign opposite to that of b */
            c = a;  fc = fa;
        }
        
    }
    /* failed! */
    *Tol = std::fabs(c-b);
    *Maxit = -1;

    return b;
}

double retrieve_lambda_bg(int k,double prob)
{
    int Maxit = 1000;
    double tol = 2.220446e-16;
    
    int lower_bound = 0;
    int upper_bound = LAMBDA_BG_UPPER_BOUND;
    
    double lower_pval = 1.0 - prob;
    
    double fa =  ppois_inv_lambda(lower_bound, k, lower_pval);
    double fb =  ppois_inv_lambda(upper_bound, k, lower_pval);
    
    
    return zeroin2(k,lower_pval,lower_bound,upper_bound,fa,fb,&tol,&Maxit);
    
    
}

