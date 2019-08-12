//
//  IntervalTree.cpp
//  BAMQC-0.5
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include "IntervalTree.h"
#include <iostream>



Interval::Interval(int g, int e, int st, int end, int t) {
    gene = g;
    start = st;
    stop = end;
    type = t;
    exon = e;
}

Interval::Interval(){}
Interval::~Interval() { }


bool itv_comp(Interval first, Interval second){
    return first.start < second.start ;
}

IntervalTree::IntervalTree(){}

IntervalTree::IntervalTree(std::vector<Interval> intervals, int depth, int minbucket, int extent_st , int extent_end,  int maxbucket){
    
    std::vector<Interval> lefts ;
    std::vector<Interval> rights ;
    int left_ext = 0;
    int right_ext = 0;
    int max_stop = 0;
    float center_pos = 0.0 ;
    Interval itv ;
    size_t i ;
    
    left = nullptr;
    right = nullptr;
    
    depth -= 1;
    if ((depth == 0 || intervals.size() < (size_t)minbucket) && intervals.size() < (size_t)maxbucket){
        itvlist = intervals;
        left = right = nullptr;
        center = -1;
        return ;
    }
    
    //if (extent_st == -1){
    // sorting the first time through allows it to get
    // better performance in searching later.

     //   quick_sort(intervals,0,intervals.size()-1);
    //}
    
    for(i=0;i<intervals.size(); i++){
        itv = intervals[i];
        if (itv.stop > max_stop) { max_stop = itv.stop; }
    }
    
    if(extent_st != -1){
        left_ext = extent_st;
        right_ext = extent_end;
    }
    else {
        
        left_ext = intervals[0].start;
        right_ext = max_stop;
    }
    //left, right = _extent or (intervals[0].start, max_stop)
    //center = intervals[len(intervals)/ 2].stop
    center_pos = (left_ext + right_ext) / 2.0;
    
    //std::cout << "left_ext " << left_ext << std::endl;
    //std::cout << "right_ext " << right_ext << std::endl;
    //std::cout << "center " << center_pos << std::endl;
    
    for(i=0;i<intervals.size();i++){
        
        itv = intervals[i];
        
        if (itv.stop < center_pos) { lefts.push_back(itv); }
        else {
            if (itv.start > center_pos) { rights.push_back(itv);} else { itvlist.push_back(itv) ;}
        }
    }
    if (lefts.empty()) {
        left = nullptr;
    }
    else{

        left = new IntervalTree(lefts,depth,minbucket,intervals[0].start,center_pos,maxbucket);
    }
    if (rights.empty()) {
        right = nullptr;
    }
    else {
        right  = new IntervalTree(rights, depth, minbucket, center_pos,right_ext,maxbucket);
    }

    center = center_pos;
    
}

void IntervalTree::clear(IntervalTree * node){
   // std::cout << "clear node" << std::endl;
    //if (node != nullptr) {
      //  std::cout << node->center << std::endl;
    //}
    if (node != nullptr) {
        clear(node->left);
        clear(node->right);
       // delete node;
    }
}

IntervalTree::~IntervalTree() {
    
    clear(left);
    clear(right);
    
/*    if (left != nullptr){
        
       if( left->left == nullptr && left->right == nullptr) {
        delete left;
       }
    
       else {
        //if (left->left != nullptr) {
            (*left).~IntervalTree();
       // }
        //if (left->right != nullptr) {
         //   (*left->right).~IntervalTree();
        //}
       }
    }
    if (right != nullptr ) {
        
       if(right->left == nullptr && right->right == nullptr) {
        delete right;
       }
    else{
       // if (right->left != nullptr) {
        //    (*right->left).~IntervalTree();
        //}
        //if (right->right != nullptr) {
         //   (*right->right).~IntervalTree();
        //}
        (*right).~IntervalTree();
     }
    }
    delete this;*/
}



std::vector<Interval> IntervalTree::find(int start, int stop){
//"""find all elements between (or overlapping) start and stop"""
    std::vector<Interval> overlapping ;
    std::vector<Interval> temp ;
    size_t i;

    if (!itvlist.empty() && stop >= itvlist[0].start) {
        for (i=0;i< itvlist.size(); i++) {
            if (itvlist[i].stop >= start && itvlist[i].start <= stop) {
                overlapping.push_back(itvlist[i]);
            }
            
        }
    }


    if( left != nullptr && start <= center) {
        temp = (*left).find(start, stop);
        overlapping.insert(overlapping.end(), temp.begin(),temp.end());
    }

    if(right != nullptr && stop >= center){
        temp = (*right).find(start, stop);

        overlapping.insert(overlapping.end(),temp.begin(),temp.end());
    }

    return overlapping ;
}

std::vector<int> IntervalTree::find_gene(int start, int stop){
    
    std::vector<int> overlapping ;
    std::vector<int> temp ;
    size_t i;
    
    if (!itvlist.empty() && stop >= itvlist[0].start) {
        for (i=0;i< itvlist.size(); i++) {
            if (itvlist[i].stop >= start && itvlist[i].start <= stop) {
                overlapping.push_back(itvlist[i].gene);
            }
            
        }
    }
    
    
    if( left != nullptr && start <= center) {
        temp = (*left).find_gene(start, stop);
        overlapping.insert(overlapping.end(), temp.begin(),temp.end());
    }
    
    if(right != nullptr && stop >= center){
        temp = (*right).find_gene(start, stop);
        
        overlapping.insert(overlapping.end(),temp.begin(),temp.end());
    }
    
    return overlapping ;
}


/*
int main(){
    std::vector<Interval> ll;
    IntervalTree *Idx;
    std::vector<Interval> res;
    std::vector<std::string> res_gene;
    int i;
    
    
    ll.push_back(Interval("gid1",100,300,"cds"));
    ll.push_back(Interval("gid2",10,100,"cds"));
    ll.push_back(Interval("gid3",1000,3000,"cds"));
    ll.push_back(Interval("gid4",500,800,"cds"));
    ll.push_back(Interval("gid5",200,700,"cds"));
    
    std::cout << "before create idx" << std::endl;
    
    Idx = new IntervalTree(ll,16, 2, -1, -1,  3);
    
    std::cout << "after create idx" << std::endl;
    
    res = (*Idx).find(20,600);
    res_gene = (*Idx).find_gene(20,600);
    
    for (i = 0; i<res.size(); i++) {
        std::cout << res[i].gene << std::endl;
    }
    
    for (i = 0; i<res_gene.size(); i++) {
        std::cout << res[i].gene << std::endl;
    }
}
*/
