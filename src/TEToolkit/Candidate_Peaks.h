//
//  Candidate_Peaks.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/22/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____Candidate_Peaks__
#define __TEToolkit_c____Candidate_Peaks__

#include <stdio.h>
//#include "IntervalTree.h"
//#include "Constants.h"
#include "GeneFeatures.h"

class Candidate_Peaks{
    
public:
    
    Candidate_Peaks(std::string Peakfilename);
    ~Candidate_Peaks();
    
    double get_count(int g);
    
    double get_pval(int g);
    double get_fe(int g);
    int get_length(int g);
    int get_start(int g);
    int get_end(int g);
    std::string get_chrom(int g);
    
    
    int get_numofpeaks();
    int get_ovp_peaks(std::string chrom,std::vector<std::pair<int,int> > itv_list, int shiftsize,std::string strand);
    

    std::vector<int> peak_length;
    
private:
    std::vector<double> peak_reads;
    std::vector<int> peak_start;
    std::vector<int> peak_end;
    std::vector<std::string> peak_chrom;
    
    
    std::vector<double> peak_pval;
    std::vector<double> peak_fe;
    std::vector<int> peak_summit;
    
    chrom_itvTree_Dict idx_peaks ;
    
    void read_peaks(std::string peak_filename) ;
    void build_tree(std::map<std::string, std::vector<std::pair<int,std::string> > > peak_chrom_maps);
    
};



#endif /* defined(__TEToolkit_c____Candidate_Peaks__) */
