//
//  Pileup.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 6/7/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//
//Modified from MACS2

#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "Pileup.h"
#include "myLog.h"

// general function to calculate maximum between two arrays.
std::pair<std::vector<int>,std::vector<double> > max_over_two_pv_array ( std::pair<std::vector<int>,std::vector<double> > tmparray1, std::pair<std::vector<int>,std::vector<double> > tmparray2 )
{
    /*"""Merge two position-value arrays. For intersection regions, take
    the maximum value within region.

    """*/
    
    //int pre_p; //, p1, p2;
    //float v1, v2;
    std::vector<int> a1_pos, a2_pos, ret_pos;
    std::vector<double> a1_v, a2_v, ret_v;
    
    long l1, l2,  i1, i2;

    a1_pos = tmparray1.first;
    a1_v = tmparray1.second;
    
    a2_pos = tmparray2.first;
    a2_v = tmparray2.second;
    

    l1 = a1_pos.size();
    l2 = a2_pos.size();

    i1 = 0;
    i2 = 0;
    //I = 0;

   // pre_p = 0;           // remember the previous position in the new bedGraphTrackI object ret

    while (i1 < l1 and i2 < l2)
    {
        if (a1_pos[i1] < a2_pos[i2]){
            // clip a region from pre_p to p1, then set pre_p as p1.
            ret_pos.push_back(a1_pos[i1]);
            ret_v.push_back(std::max( a1_v[i1], a2_v[i2] ));

          //ADD  pre_p = a1_pos[i1];
            // call for the next p1 and v1
            i1 += 1;
        }
        else {
        if (a1_pos[i1] > a2_pos[i2])
        {
            // clip a region from pre_p to p2, then set pre_p as p2.
            ret_pos.push_back(a2_pos[i2]);
            ret_v.push_back(std::max( a1_v[i1], a2_v[i2] ));

            //I += 1;
         //ADD   pre_p = a2_pos[i2];
            // call for the next p1 and v1
            i2 += 1;
        }
        else
        {
            // from pre_p to p1 or p2, then set pre_p as p1 or p2.
            ret_pos.push_back(a1_pos[i1]);
            ret_v.push_back(std::max( a1_v[i1], a2_v[i2]));

            //I += 1;
         //ADD   pre_p = a1_pos[i1];
            // call for the next p1, v1, p2, v2.
            i1 += 1;
            i2 += 1;
        }
        }
    }
    
    return std::pair<std::vector<int>,std::vector<double> >(ret_pos, ret_v);
    
}



void fix_coordinates(std::vector<int> &poss, int rlength)
{
    if( poss.size() == 0) 
    {
       return ;
    }
    for(size_t i=0; i < poss.size() ; i++)
    {
        if (poss[i] < 0){
            poss[i] = 0;
        }
        else{
            break;
        }
    }
    for(size_t i = poss.size() -1; i >= 0 ; i--)
    {
        if (poss[i] > rlength){
            poss[i] = rlength;
        }
        else{
            break;
        }
    }
}

std::pair<std::vector<int>, std::vector<double> > se_all_in_one_pileup ( std::vector<int> plus_tags, std::vector<int> minus_tags, int five_shift, int three_shift, int rlength, double scale_factor, double baseline_value )
{
    /*"""Return pileup given 5' end of fragment at plus or minus strand
    separately, and given shift at both direction to recover a
    fragment. This function is for single end sequencing library
    only. Please directly use 'quick_pileup' function for Pair-end
    library.
    
    It contains a super-fast and simple algorithm proposed by Jie
    Wang. It will take sorted start positions and end positions, then
    compute pileup values.

    It will return a pileup result in similar structure as
    bedGraph. There are two python arrays:
    
    [end positions, values] or '[p,v] array' in other description for
    functions within MACS2.

    Two arrays have the same length and can be matched by index. End
    position at index x (p[x]) record continuous value of v[x] from
    p[x-1] to p[x].

    """ */
    
    long i_s, i_e,  I;
    int  p, pre_p, pileup;
    
    std::vector<int> start_poss, end_poss;
    
    //pointers are used for numpy arrays
   // int * start_poss_ptr;
   // int * end_poss_ptr;
   // int * ret_p_ptr;     // pointer for position array
    //float * ret_v_ptr ;    // pointer for value array
        

    for (size_t i= 0; i < plus_tags.size(); i++) {
        start_poss.push_back(plus_tags[i] - five_shift);
        end_poss.push_back(plus_tags[i] + three_shift);
    }
    for (size_t i=0; i < minus_tags.size(); i++) {
        start_poss.push_back(minus_tags[i] - three_shift);
        end_poss.push_back(minus_tags[i] + five_shift);
    }
     debug( "start se_all_in_one_pileup " + std::to_string(start_poss.size()) + "\t" + std::to_string(rlength) );
    // sort
    std::sort(start_poss.begin(),start_poss.end());
    std::sort(end_poss.begin(),end_poss.end());
    
    // fix negative coordinations and those extends over end of chromosomes
    fix_coordinates(start_poss, rlength);
    fix_coordinates(end_poss, rlength);
    
    //lx = start_poss.size();
   // start_poss_ptr = <int32_t *> start_poss.data
   // end_poss_ptr = <int32_t *> end_poss.data

    std::vector<int> ret_p; //(2*lx,0);
    std::vector<double> ret_v; //( 2 * lx, 0.0 );

    int lx = start_poss.size();
   // ret_p_ptr = <int32_t *> ret_p.data
   // ret_v_ptr = <float32_t *> ret_v.data

    //std::pair<std::vector<int>,std::vector<double> >tmp (ret_p,ret_v);//for (endpos,value)
    
    i_s = 0 ;                        // index of start_poss
    i_e = 0 ;                        // index of end_poss
    I = 0 ;

    pileup = 0;
    
    if (start_poss.size() == 0)
    {
        return  std::pair<std::vector<int>,std::vector<double> >(ret_p,ret_v);;
    }
    
    pre_p = std::min(start_poss[0],end_poss[0]);

    if (pre_p != 0){
        // the first chunk of 0
        ret_p.push_back(pre_p);
        ret_v.push_back(std::max(0.0,baseline_value)) ;
        I += 1;
    }
    //pre_v = pileup;

    if (start_poss.size() != end_poss.size()) {
        error("error in pileup, different number of start_pos and end_pos!");
        std::exit(1);
    }

   
    while (i_s < lx and i_e < lx){
        if (start_poss[i_s] < end_poss[i_e])
        {
            p = start_poss[i_s];
            if ( p != pre_p){
                ret_p.push_back(p) ;
                ret_v.push_back(std::max(pileup * scale_factor, baseline_value));
                
                pre_p = p;
            }
            pileup += 1;
            i_s += 1;
            //start_poss_ptr += 1
        }
        else {
        if (start_poss[i_s] > end_poss[i_e])
        {
            p = end_poss[i_e];
            if (p != pre_p){
                ret_p.push_back(p);
                ret_v.push_back(std::max(pileup * scale_factor, baseline_value));

                
                pre_p = p;
            }
            pileup -= 1;
            i_e += 1;
            //end_poss_ptr += 1;
        }
        else {
            i_s += 1;
            i_e += 1;
        }
        }
    }
    
    //debug("after while loop");
    if (i_e < lx)
    {
        // add rest of end positions
        for(int i = i_e; i < lx ; i ++ )
        {
            p = end_poss[i];
            if (p != pre_p)
            {
                ret_p.push_back(p);
                ret_v.push_back(std::max(pileup * scale_factor, baseline_value));

                I += 1;
                pre_p = p;
            }
            pileup -= 1;
            //end_p += 1;
        }
    }

    return std::pair<std::vector<int>,std::vector<double> >(ret_p,ret_v);
}


std::pair<std::vector<int>,std::vector<double> > quick_pileup ( std::vector<pos_t> poss, double scale_factor, double baseline_value )
{
    /*"""Return pileup given plus strand and minus strand positions of fragments.
    
    A super-fast and simple algorithm proposed by Jie Wang. It will
    take sorted start positions and end positions, then compute pileup
    values. 

    It will return a pileup result in similar structure as
    bedGraph.

    Two arrays have the same length and can be matched by index. End
    position at index x (p[x]) record continuous value of v[x] from
    p[x-1] to p[x].

    """*/
    
    
    int i_s, i_e,  I;
    int p, pre_p, pileup;
    
    int ls = poss.size();
    //int le = end_poss.size();
    //int l = ls + le;
    std::vector<int> ret_p;
    std::vector<double> ret_v;
    
    std::pair<std::vector<int>,std::vector<double> > tmp;
    // pointers are used for numpy arrays

    i_s = 0;
    i_e = 0;
    I = 0;

    pileup = 0;
    
    if (ls == 0){
        return tmp;
    }
    pre_p = std::min(poss[0].start, poss[0].end);

    if (pre_p != 0)
    {
        ret_p.push_back(pre_p);
        ret_v.push_back(std::max(0.0,baseline_value)) ;

        I += 1;
        i_s +=1;
        i_e += 1;
    }
    
    //pre_v = pileup;
    
    while( i_s < ls and i_e < ls)
    {
        if (poss[i_s].start < poss[i_e].end)
        {
            p = poss[i_s].start;
            if (p != pre_p)
            {
                ret_p.push_back(p);
                ret_v.push_back(std::max(pileup * scale_factor, baseline_value));
                
                i_s += 1;
                I += 1;
                pre_p = p;
            }
            pileup += 1;

            i_s += 1;
        }
        if (poss[i_s].start > poss[i_e].end)
        {
            p = poss[i_e].end;
            if (p != pre_p)
            {
                ret_p.push_back(p);
                ret_v.push_back(std::max(pileup * scale_factor, baseline_value));
                I += 1;
                pre_p = p;
            }
            pileup -= 1;
            i_e += 1;
        }
        if (poss[i_s].start == poss[i_e].end){
            i_s += 1;
            i_e += 1;
        }
    }

    if (i_e < ls)
    {
        // add rest of end positions
        for (int i = i_e; i< ls; i++)
        {
            p = poss[i].end;
            //for p in minus_tags[i_e:]:
            if (p != pre_p)
            {
                ret_p.push_back(p);
                ret_v.push_back(std::max(pileup * scale_factor, baseline_value));
                I += 1;
                pre_p = p;
            }
            pileup -= 1;
        }
    }
    tmp.first = ret_p;
    tmp.second = ret_v;
    return tmp;
}
