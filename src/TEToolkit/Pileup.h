//
//  Pileup.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 6/7/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____Pileup__
#define __TEToolkit_c____Pileup__

#include <stdio.h>
#include "ShortRead.h"

std::pair<std::vector<int>, std::vector<double> > se_all_in_one_pileup ( std::vector<int> plus_tags, std::vector<int> minus_tags, int five_shift, int three_shift, int rlength, double scale_factor, double baseline_value );

std::pair<std::vector<int>,std::vector<double> > max_over_two_pv_array ( std::pair<std::vector<int>,std::vector<double> > tmparray1, std::pair<std::vector<int>,std::vector<double> > tmparray2 );

std::pair<std::vector<int>,std::vector<double> > quick_pileup ( std::vector<pos_t> poss, double scale_factor, double baseline_value );

#endif /* defined(__TEToolkit_c____Pileup__) */
