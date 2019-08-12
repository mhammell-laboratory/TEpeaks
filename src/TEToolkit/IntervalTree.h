//
//  IntervalTree.h
//  BAMQC-0.5
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef __BAMQC_0_5__IntervalTree__
#define __BAMQC_0_5__IntervalTree__

#include <stdio.h>
#include <string>
#include <vector>

#include "Constants.h"



class Interval{
public:
    int  start;
    int stop;
    //std::string gene;
    int gene;
    //std::string type;
    int type;
    int exon;

    Interval();
    Interval(int g, int exon,int st, int end, int t);
    ~Interval();
};





class IntervalTree{
public :
    std::vector<Interval> itvlist;
    IntervalTree* left ;
    IntervalTree* right ;
    float center ;

    IntervalTree();
    IntervalTree(std::vector<Interval> intervals, int depth = DEPTH, int minbucket= MIN_BUCKET, int extent_st=-1, int extent_end = -1 , int maxbucket= MAX_BUCKET);
    ~IntervalTree();
    std::vector<Interval> find(int start, int stop);
    std::vector<int>  find_gene(int start, int stop);
private:
    void clear(IntervalTree * node);

};

extern "C" {
    bool itv_comp(Interval first, Interval second);
}


#endif /* defined(__BAMQC_0_5__IntervalTree__) */
