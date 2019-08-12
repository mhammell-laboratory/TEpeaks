//
//  PeakDetect.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 5/3/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____PeakDetect__
#define __TEToolkit_c____PeakDetect__

#include <stdio.h>
#include "Parser.h"
#include "ShortRead.h"
#include "PeakIO.h"


/* A PeakInfo type is a: dictionary
key value: chromosome

items: PeakContent (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment) <-- tuple type
*/


//typedef std::map<int,std::vector<PeakContent> > PeakInfo;

class PeakDetect
{
public:
    opt_t * opt;
    
    ShortRead * treat ;
    ShortRead * control ;
    
    double ratio_treat2control ;
    PeakIO * peaks ;
    
    //bool PE_MODE ;
    //std::string scoretrack;
    
    
  //  double log_pvalue = opt.log_pvalue;    // -log10pvalue
  //  double log_qvalue = opt.log_qvalue;    // -log10qvalue
    
    //PeakInfo peakinfos; //?
    int d ;
    //int gsize ;
    //int end_shift ;
    //bool nolambda ;
    //int sregion ;
    //int lregion ;
    
    PeakDetect(ShortRead * treat,ShortRead *control,opt_t *options);
    ~PeakDetect();
    
    void call_peaks();
    void candidate_peakregions();
    //void filter_fc(double fc_low);
    
private:
 //ADD   void __call_peaks_wo_control ();
    void __call_peaks_w_control (bool uniqOnly);
};
#endif /* defined(__TEToolkit_c____PeakDetect__) */
