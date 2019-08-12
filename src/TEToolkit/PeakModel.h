//
//  PeakModel.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 4/29/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____PeakModel__
#define __TEToolkit_c____PeakModel__

#include <stdio.h>
#include <string>
//#include <vector>
//#include <map>

#include "ShortRead.h"
#include "Parser.h"

//typedef std::map<std::string, std::vector<std::string> > paired_peaks_Dict ;
typedef std::map<int,std::vector<int> > PairedPeaks_Dict;

class PeakModel {
    
public :
    ShortRead *treatment;
    double gz;
    int max_pairnum;
    int umfold;
    int lmfold;
    int bw;
    int tag_expansion_size;
    //object info, debug, warn, error;
    std::string summary;
    
    //public np.ndarray plus_line, minus_line, shifted_line;
    int d ;
    int scan_window ;
    int min_tags;
    int max_tags ;
    int peaksize ;
    std::vector<int> alternative_d ;
    bool NotEnoughPairs;
    std::vector<int> xcorr;
    std::vector<double> ycorr;
    std::vector<int> plus_line;
    std::vector<int> minus_line;
    
    
    PeakModel(ShortRead *treatment, opt_t options, int max_pairnum);
    bool build() ;
    ~PeakModel();
    
private:
    void __paired_peak_model( PairedPeaks_Dict paired_peakpos);
    void __model_add_line ( std::vector<int> pos1, std::vector<int> tags,  std::vector<int>& start, std::vector<int>& end);
    void __count ( std::vector<int> start, std::vector<int> end, std::vector<int>& line );
    
     void __paired_peaks (PairedPeaks_Dict &paired_peakpos);
    
    std::vector<int> __find_pair_center (std::vector<std::pair<int,int> > pluspeaks, std::vector<std::pair<int,int> > minuspeaks);
    
    std::vector<std::pair<int, int> > __naive_find_peaks (std::vector<int> taglist, int size, int plus_strand=1 );
    
    int __naive_peak_pos (std::vector<int> pos_list, int plus_strand=1 );
    

};
#endif /* defined(__TEToolkit_c____PeakModel__) */
