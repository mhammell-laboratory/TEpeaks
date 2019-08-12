//
//  EM_TEpeaks_BAM.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/22/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#include "EM_TEpeaks.h"

//#include <malloc.h>
#include <string>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <ctime>
#include <sstream>
#include <tuple>

#include "htslib/sam.h"
//#include "IntervalTree.h"
#include "Candidate_Peaks.h"
#include "CallerFromAlignments.h"
#include "cStatistics.h"

#include "EMestimate_reads.h"
#include "zeroin.h"
#include "myLog.h"
//#include "Parser.h"


void prepare_EM_structure(Candidate_Peaks *peakIdx,std::vector<double> peak_reads, MULTI_READ_MAP * multi_read_mapTo_pids, std::vector<int> & multiAlgn_To_multiRead,std::vector<double> & multiAlgn_weight, std::vector<int> & multiAlgn_To_peakdID,std::vector<int> & multiAlgn_start)
{
    //std::map<int,std::vector<int> >::iterator 
    MULTI_READ_MAP::iterator itr;
    //multi-aligment -> peakID
    int multiAlgn_id = 0;
    
    std::cout << "multi read map to peak " << multi_read_mapTo_pids->size() << std::endl;
    
    for (itr = multi_read_mapTo_pids->begin(); itr != multi_read_mapTo_pids->end(); itr ++)
    {
        //int readID = itr->first;
        //std::cout << "multi read " << readID << std::endl;
        std::vector<int> pidlist = itr->second.pidlist;
        std::vector<int> startlist = itr->second.startlist;
        
        double total = 0.0; // std::accumulate(counts_.begin(), counts_.end(), 0.0);;
        for (size_t i=0; i < pidlist.size(); i++) {
            //std::cout << "pid " << pidlist[i] << std::endl;
            //std::cout << "start " << startlist[i] << std::endl;
            total += 1.0 * peak_reads[pidlist[i]]/peakIdx->get_length(pidlist[i]);
            
          //  std::cout <<  "#1 align_id = " << multiAlgn_id << " pid = " << pidlist[i] << std::endl;
            
            multiAlgn_id += 1;
            multiAlgn_To_peakdID.push_back(pidlist[i]);
            multiAlgn_start.push_back(startlist[i]);
          //  std::cout <<  "#1 align_id = " << multiAlgn_id << " pos = " << startlist[i] << std::endl;
            
        }
        multiAlgn_To_multiRead.push_back(multiAlgn_id);
        for (size_t i=0; i < pidlist.size(); i++) {
            double w = 1.0 *peak_reads[pidlist[i]]/(peakIdx->get_length(pidlist[i])*total);
            multiAlgn_weight.push_back(w);
        }
    }
  //  if (multiAlgn_id > 0) {
  //      multiAlgn_To_multiRead.push_back(multiAlgn_id);
   // }
}

//equaly weight multi-aignments
//the input peak_reads contains unique reads only
void EqualWeight_reads(std::vector<int> multiAlgn_To_multiRead, std::vector<double> peak_reads, std::vector<double> & peak_reads_Prime,std::vector<int> multiAlgn_To_peakID,std::vector<double> & multiAlgn_weight)
{
    //debug("in EqualWeight_reads");
    int multiAlgn_first_idx = 0;
    //std::vector<double> tmp_peak_reads(peak_reads.size(),0.0);
    for (size_t i= 0; i < peak_reads.size(); i++) {
        peak_reads_Prime[i] = peak_reads[i];
    }
    
//multi-alignment of the same read are saved consequtively
    for (size_t i=0; i< multiAlgn_To_multiRead.size(); i++) {
        //loop through all multi-alignments of each multi-read
        double w = 1.0 / (multiAlgn_To_multiRead[i] - multiAlgn_first_idx);
        
        for (int j= multiAlgn_first_idx; j < multiAlgn_To_multiRead[i]; j++) {
            int pid = multiAlgn_To_peakID[j];
            
            peak_reads_Prime[pid] += w;
            multiAlgn_weight.push_back(w);
        }
        multiAlgn_first_idx = multiAlgn_To_multiRead[i];
    }
}


void pair_treat_ctrl ( std::pair<std::vector<int>,std::vector<double> > treat_pv, std::pair<std::vector<int>,std::vector<double> > ctrl_pv, pos_tc_t &pos_treat_ctrl_pair )
{
    /*"""*private* Pair treat and ctrl pileup for each peak.
     treat_pv and ctrl_pv are [np.ndarray, np.ndarray].
     return [p, t, c] list, each element is a vector.
     """*/
    
    int lt = treat_pv.first.size();
    int lc = ctrl_pv.first.size();
    
    //std::vector<pos_tc_t> ret;
    
    //int pre_p = 0;
    //int index_ret = 0;
    int it = 0;
    int ic = 0;
    
    while (it < lt and ic < lc){
        if (treat_pv.first[it] < ctrl_pv.first[ic])
        {
            // clip a region from pre_p to p1, then set pre_p as p1.
            pos_treat_ctrl_pair.pos.push_back(treat_pv.first[it]);
            pos_treat_ctrl_pair.treat_v.push_back(treat_pv.second[it]);
            pos_treat_ctrl_pair.ctrl_v.push_back(ctrl_pv.second[ic]);
            
            //pre_p = treat_pv.first[it];
            //index_ret += 1;
            // call for the next p1 and v1
            it += 1;
        }
        else {
            if(treat_pv.first[it] > ctrl_pv.first[ic]){
                // clip a region from pre_p to p2, then set pre_p as p2.
                
                pos_treat_ctrl_pair.pos.push_back(ctrl_pv.first[ic]);
                pos_treat_ctrl_pair.treat_v.push_back(treat_pv.second[it]);
                pos_treat_ctrl_pair.ctrl_v.push_back(ctrl_pv.second[ic]);
                
                //pre_p = ctrl_pv.first[ic];
                // index_ret += 1;
                // call for the next p2 and v2
                ic += 1;
            }
            else{
                // from pre_p to p1 or p2, then set pre_p as p1 or p2.
                pos_treat_ctrl_pair.pos.push_back(treat_pv.first[it]);
                pos_treat_ctrl_pair.treat_v.push_back(treat_pv.second[it]);
                pos_treat_ctrl_pair.ctrl_v.push_back(ctrl_pv.second[ic]);
                
                // pre_p = treat_pv.first[it];
                //index_ret += 1;
                // call for the next p1, v1, p2, v2.
                it += 1;
                ic += 1;
                
            }
        }
    }
    //return ret;
}

//code position and weight into one double value
std::pair<int,double> get_pos_id(double pos_w,std::vector<double> multi_weights,std::vector<int> maid_list)
{
    //get the first start position
    int pos = (int) std::floor(pos_w);
    double pre_v = pos_w - pos;
    //int id = -1;
    double w = 1.0;
    
 /*   if(pre_v != 0)//multi_read
    {
        id = (int) (1.0/pre_v) - 2 ;
        w = multi_weights[id];
    }
    else { //uniq_read
        w = 1.0;
    }*/
    if(pre_v != 0)//multi_read
    {
        int idx = (int) (1.0/pre_v) - 2 ;
        //info("idx = " + std::to_string(idx));
        int ma_id = maid_list[idx];
        if ( ma_id > multi_weights.size())
        {
            
            info("ma_id = " + std::to_string(ma_id));
            info("multi_weights.size = " + std::to_string(multi_weights.size()));
            info("id is not in multi_weights in get_pos_id !!");
            //std::exit(1);
            w = 0.0;
        }
        if(ma_id < 0)
        {
            info("id is less than 0!!");
            info(std::to_string(pre_v));
            info("pos_w " + std::to_string(pos_w));
            w = 0.0;
            //std::exit(1);
        }
        w = multi_weights[ma_id];
    }
    else { //uniq_read
        w = 1.0;
    }
    
    return std::pair<int,double> (pos,w);
    
}

void pileup_peak ( std::vector<std::tuple<int,double> > poss, int fragsize, std::pair<std::vector<int>,std::vector<double> > &ret )
{
    /*"""Return pileup given start positions of fragments of both unique reads and multi-reads.
     
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
    int p, pre_p;
    double pileup;
    
    int ls = poss.size();
    //int le = end_poss.size();
    //int l = ls + le;
   // std::vector<int> ret_p;
   // std::vector<double> ret_v;
    
    //info("poss size " + std::to_string(ls));
    //std::pair<std::vector<int>,std::vector<double> > tmp;
    // pointers are used for numpy arrays
    
    i_s = 0;
    i_e = 0;
    I = 0;
    
    pileup = 0;
    
    if (ls == 0){
        return ;
    }
    
    //get the first start position
    //std::pair<int,double> pos_id = get_pos_id(poss[0],multi_weights,maid_list);
    std::tuple<int,double> pos_id = poss[0];
    //info("after get pos_id 0");
    pre_p = std::get<0>(pos_id);
    // pre_w = pos_id.second;
    
    
    ret.first.push_back(pre_p );
    ret.second.push_back(0.0);
    
    I += 1;
    
    std::vector<std::tuple<int,double> > pos_end;
    
    for (auto k : poss ) {
        
        int k_pos = std::get<0>(k) + fragsize;
        double k_w = std::get<1>(k);
        
        pos_end.push_back(std::tuple<int,double>(k_pos,k_w));
        
    }
    
    //pre_v = pileup;
    
    while( i_s < ls and i_e < ls)
    {
       
        std::tuple<int,double> pos_id_s = poss[i_s];
        
        std::tuple<int,double> pos_id_e = pos_end[i_e];
        
        int p_start = std::get<0>(pos_id_s);
        int p_end = std::get<0>(pos_id_e);
        
        if (p_start < p_end)
        {
            //p = pos_id_s.first;
            p = std::get<0>(pos_id_s);
            
            //double w = pos_id_s.second;
            double w = std::get<1>(pos_id_s);
            
            if (p != pre_p)
            {
                ret.first.push_back(p);
                ret.second.push_back(pileup );
                
                I += 1;
                pre_p = p;
            }
            pileup += w;
            
            i_s += 1;
        }
        
        if (p_start > p_end)
        {
            
            p = std::get<0>(pos_id_e);
            double w = std::get<1>(pos_id_e);
            
            if (p != pre_p)
            {
                ret.first.push_back(p);
                ret.second.push_back(pileup);
                I += 1;
                pre_p = p;
            }
            pileup -= w;
            
            i_e += 1;
        }
        if (p_start == p_end){
            
            double pre_w = std::get<1>(pos_id_e);
            
            double w = std::get<1>(pos_id_s);
            p = std::get<0>(pos_id_s);
            
            if(p != pre_p)
            {
                ret.first.push_back(p);
                ret.second.push_back(pileup);
            }
            pileup -= pre_w;
            
            pileup += w;
            
            i_s += 1;
            i_e += 1;
        }
    }
   
    if (i_e < ls)
    {
        // add rest of end positions
        for (int i = i_e; i< ls; i++)
        {
            
            std::tuple<int,double> pos_id = pos_end[i];
            
            p = std::get<0>(pos_id);
            double w = std::get<1>(pos_id);
            
            if (p != pre_p)
            {
                ret.first.push_back(p);
                ret.second.push_back(pileup);
                I += 1;
                pre_p = p;
            }
            pileup -= w;
        }
    }
    //tmp.first = ret_p;
    //tmp.second = ret_v;
   // return tmp;
}

//void write_multi_read_bedGraph_for_a_chromosome(std::string chrom,pos_tc_t pos_treat_ctrl_multi_pair, std::ofstream &bedGraph_treat_f, std::ofstream  &bedGraph_ctrl_f)
void write_multi_read_bedGraph_for_a_chromosome(std::string chrom,pos_tc_t pos_treat_ctrl_multi_pair, std::string outdir, std::string prj_name)
{
    //"""Write treat/control values for a certain chromosome into a
    //    specified file handler.
    //    """
    std::vector<int> pos_array;
    std::vector<double> treat_array;
    std::vector<double> ctrl_array;
    int l,lt,lc;
    int p,pre_p_t,pre_p_c; //current position, previous position for treat, previous position for control
    double pre_v_t, pre_v_c, v_t, v_c; // previous value for treat, for control, current value for treat and control
    
    //double denominator = 1; //
    
    pos_array = pos_treat_ctrl_multi_pair.pos;
    treat_array = pos_treat_ctrl_multi_pair.treat_v;
    ctrl_array = pos_treat_ctrl_multi_pair.ctrl_v;
    
    l = pos_array.size();
    lt = treat_array.size();
    lc = ctrl_array.size();
    if(l != lt || l != lc)
    {
        info("arrays have different sizes !!");
        std::cout << l << std::endl;
        std::cout << lt << std::endl;
        std::cout << lc << std::endl;
        
        std::exit(1);
    }
    if(l ==0)
    {
        return;
    }
    std::ofstream  bedGraph_treat_f;
    std::ofstream  bedGraph_ctrl_f;
    std::string  bedGraph_filename_prefix = outdir+"/" + prj_name;
    std::string bedGraph_treat_filename = bedGraph_filename_prefix + "_TREAT_multi.bdg";
    std::string bedGraph_control_filename = bedGraph_filename_prefix + "_CTRL_multi.bdg";
    
    bedGraph_treat_f.open(bedGraph_treat_filename,std::ofstream::app);
    bedGraph_ctrl_f.open( bedGraph_control_filename,std::ofstream::app);
    
    
    pre_p_t = pos_array[0];
    pre_p_c = pos_array[0];
    
    pre_v_t = treat_array[0];
    pre_v_c = ctrl_array[0];
    
    for(int i =1; i < l ; i ++)
    {
        v_t = treat_array[i];
        v_c = ctrl_array[i];
        p = pos_array[i];
        if (std::abs(pre_v_t - v_t) > 1e-5) // precision is 5 digits
        {
            if(pre_v_t > 1e-5)
            {
                //pre_v_t = 1e-5;
            
                bedGraph_treat_f << chrom << "\t" << pre_p_t << "\t" << p << "\t" << pre_v_t << "\n";
            }
        
            pre_v_t = v_t;
            pre_p_t = p;
        }
        if ( std::abs(pre_v_c - v_c) > 1e-5) // precision is 5 digits
        {
            if(pre_v_c > 1e-5)
            {
                //pre_v_c = 1e-5;
            
            
            bedGraph_ctrl_f << chrom << "\t" << pre_p_c << "\t" << p << "\t" << pre_v_c << "\n";
            }
            pre_v_c = v_c;
            pre_p_c = p;
        }
        
    }
    p = pos_array[l-1];
    if(pre_v_t > 1e-5)
    {
    bedGraph_treat_f << chrom << "\t" << pre_p_t << "\t" << p << "\t" << pre_v_t << "\n";
    }
    if(pre_v_c > 1e-5)
    {
    bedGraph_ctrl_f << chrom << "\t" << pre_p_c << "\t" << p << "\t" << pre_v_c << "\n";
    }
    
    bedGraph_treat_f.close();
    bedGraph_ctrl_f.close();
}

/*double roundTop(double num, int precision)
{
    
    return floor(num * pow(10.0f,precision) + .5f)/pow(10.0f,precision);
}*/

void filter_peaks(std::map<int, peak_align_t> peak_multiAlign_map,
                  std::map<int, std::vector<int> > * treat_peak_uniqReads_startPos,
                  std::map<int, std::vector<int> > * ctrl_peak_uniqReads_startPos,
                  std::vector<int> peaks_with_multiReads,
                  
                  Candidate_Peaks * peakIdx,
                  std::vector<double> IP_peakReads,
                  std::vector<double> Input_peakReads,
                  
                  std::vector<double> treat_multiAlgn_weight,
                  std::vector<int> treat_multiAlgn_start,
                  std::vector<int> ctrl_multiAlgn_start,
                  std::vector<double> ctrl_multiAlgn_weight,
                  opt_t options)
{
    
    try {
        //int fragsize = options.fragsize;
        std::string project_name = options.project_name;
       // double pvalCutoff = std::pow(10,options.log_pvalue);
        std::string data_outdir = options.data_outdir;
        
        //ADD int window_size = fragsize * 2;
        std::ofstream ROUT;
        
        
        ROUT.open (data_outdir +"/" +project_name+ "_multiRead_peaks.txt", std::ofstream::out);
        
        ROUT << "chrom\t" << "start\t" << "end\t" << "IP_tags\t" << "Input_tags\t" << "pileup\t" << "fold enrichment\t" << "p-value\n";
        
        if(IP_peakReads.size() != peaks_with_multiReads.size())
        {
            error("inconsistency of peak numbers");
            std::exit(1);
        }
         info("loop through all peaks");
        //std::ofstream f;
        //f.open("tmp",std::ofstream::out);
        std::string p_chrom = "";
        int p_start = 0;
        int p_end = 0;
        double p_summit_pscore = 0.0;
        int p_ip_tag = 0;
        int p_input_tag = 0;
        int p_summit_pileup = 0;
        double p_fe = 0.0;

        
        for (size_t i =0; i < IP_peakReads.size(); i++) {
            
           // info("i = " + std::to_string(i));
            
            if (peaks_with_multiReads[i] == 0) {
                continue;
            }
            
            double ip_tag = IP_peakReads[i];
            double input_tag = Input_peakReads[i];
          //  info(" peak start " + std::to_string(peakIdx->get_start(i)));
           // info(" peak end " + std::to_string(peakIdx->get_end(i)));
            //info("ip_tag = " + std::to_string(ip_tag * options.treat_scale));
            //info("input_tag = " + std::to_string(input_tag * options.ctrl_scale));

            double fe =  (ip_tag + 0.1)* options.treat_scale/((0.1 + input_tag)* options.ctrl_scale);
            
            double pileup = 1.0 * ip_tag * options.treat_scale;
            
           // std::vector<double> treat_pos_id_array;
            std::vector<std::tuple<int,double> > treat_pos_id_array;
            //std::vector<double> ctrl_pos_id_array;
            std::vector<std::tuple<int,double> > ctrl_pos_id_array;
            
            //std::vector<double> treat_multi_pos_id_array; //for output bedgraph
            //std::vector<double> ctrl_multi_pos_id_array; //for output bedgraph
            
            std::vector<std::tuple<int, double> > treat_multi_pos_id_array;
            std::vector<std::tuple<int, double> > ctrl_multi_pos_id_array;
           // info("options.fe = " + std::to_string(options.fe));
            
            //if (pileup > options.pileup && fe > options.fe )
            if(pileup > 5 && fe > 1.5)
            {
               //  info("fe = " + std::to_string(fe));
               //  info("peak id = " +std::to_string(i));
                //copy unique reads start position of treatment sample
                if(treat_peak_uniqReads_startPos->find(i) != treat_peak_uniqReads_startPos->end())
                {
                //info("find treat uniq reads");
                    
                   // std::copy(treat_peak_uniqReads_startPos->at(i).begin(), treat_peak_uniqReads_startPos->at(i).end(), std::back_inserter(treat_pos_id_array));
                   // int nn = treat_pos_id_array.size() - 1;
                    
                    for(auto pos : treat_peak_uniqReads_startPos->at(i))
                    {
                        treat_pos_id_array.push_back(std::make_tuple(pos,1.0));
                        //info("treat_pos_id_array[i] = " + std::to_string(treat_pos_id_array[i]));
                        //info("treat_pos_id_array[size] = " + std::to_string(treat_pos_id_array[nn]));
                    }
                }
                //copy unique read start positiion of ctrl sample
                if(ctrl_peak_uniqReads_startPos->find(i) != ctrl_peak_uniqReads_startPos->end())
                {
                // info("find ctrl uniq reads");
                //std::copy(ctrl_peak_uniqReads_startPos->at(i).begin(), ctrl_peak_uniqReads_startPos->at(i).end(),std::back_inserter(ctrl_pos_id_array));
                    //info("ctrl_pos_id_array[0] = " + std::to_string(ctrl_pos_id_array[0]));
                    for(auto pos : ctrl_peak_uniqReads_startPos->at(i))
                    {
                        ctrl_pos_id_array.push_back(std::make_tuple(pos,1.0));
                    }
                }
                
                int frag_size = options.fragsize;
                std::vector<int> treat_maid_list;
                std::vector<int> ctrl_maid_list;
                
                if(peak_multiAlign_map.find(i) != peak_multiAlign_map.end())
                {
                
                 //info("find peak multiAlign map");
                std::vector<int> treat_multi = peak_multiAlign_map[i].treat_multi_align_list;
                std::vector<int> ctrl_multi  = peak_multiAlign_map[i].ctrl_multi_align_list;
                
                //loop through all multi_alignment mapped to this peak in IP sample
                    int index = 0;
                for (auto ma_id  : treat_multi ) {
                   // info("ma id = " + std::to_string(ma_id));
                    int pos = treat_multiAlgn_start[ma_id];
                   // info("pos = " + std::to_string(pos));
                    
                    double w =treat_multiAlgn_weight[ma_id];
                    if(w >= 0.01){
                        treat_maid_list.push_back(ma_id);
                        //double code = 1.0 * pos + 1.0/(index+2);
                        
                        //treat_pos_id_array.push_back(code);
                        treat_pos_id_array.push_back(std::make_tuple(pos,w));
                        //treat_multi_pos_id_array.push_back(code);
                        treat_multi_pos_id_array.push_back(std::make_tuple(pos,w));
                        index += 1;
                    }
                }
                //loop through control sample
                index = 0;
                for (auto ma_id  : ctrl_multi ) {
                    int pos = ctrl_multiAlgn_start[ma_id];
                    double w =ctrl_multiAlgn_weight[ma_id];
                    if(w < 0)
                    {
                        info("ctrl weight less than zero !");
                        std::exit(1);
                    }
                    if(w >= 0.01){
                        ctrl_maid_list.push_back(ma_id);
                      //  double code = 1.0 * pos + 1.0/(index+2);
                     //   info("ctrl pos = " + std::to_string(pos));
                     //   info("ctrl code = " + std::to_string(code));
                        ctrl_pos_id_array.push_back(std::make_tuple(pos,w));
                        ctrl_multi_pos_id_array.push_back(std::make_tuple(pos,w));
                        index += 1;
                    }
                }
                }
                std::sort(treat_pos_id_array.begin(),treat_pos_id_array.end());
                std::sort(ctrl_pos_id_array.begin(),ctrl_pos_id_array.end());
                std::sort(treat_multi_pos_id_array.begin(),treat_multi_pos_id_array.end());
                std::sort(ctrl_multi_pos_id_array.begin(),ctrl_multi_pos_id_array.end());
                
               
                //info("before pileup_peak1");
                std::pair<std::vector<int>,std::vector<double> > treat_pos_weight_pair;
                pileup_peak(treat_pos_id_array,frag_size,treat_pos_weight_pair);
                // info("after pileup_peak1");
                std::pair<std::vector<int>,std::vector<double> > ctrl_pos_weight_pair;
                pileup_peak(ctrl_pos_id_array,frag_size,ctrl_pos_weight_pair);
                // info("after pileup_peak2");
                //for output bedgraph of multi_reads
                std::pair<std::vector<int>,std::vector<double> > treat_multi_pos_weight_pair ;
                pileup_peak(treat_multi_pos_id_array,frag_size,treat_multi_pos_weight_pair);
                
                // info("after pileup_peak3");
                std::pair<std::vector<int>,std::vector<double> > ctrl_multi_pos_weight_pair;
                pileup_peak(ctrl_multi_pos_id_array,frag_size,ctrl_multi_pos_weight_pair);
                // info("after pileup_peak4");
                
                
                pos_tc_t pos_treat_ctrl_pair;
                pos_treat_ctrl_pair.pos.clear();
                pos_treat_ctrl_pair.treat_v.clear();
                pos_treat_ctrl_pair.ctrl_v.clear();
                
                pair_treat_ctrl(treat_pos_weight_pair,ctrl_pos_weight_pair,pos_treat_ctrl_pair);
                
                // info("after pair_treat_ctrl1");
                pos_tc_t pos_treat_ctrl_multi_pair;
                pos_treat_ctrl_multi_pair.pos.clear();
                pos_treat_ctrl_multi_pair.treat_v.clear();
                pos_treat_ctrl_multi_pair.ctrl_v.clear();
                
                pair_treat_ctrl(treat_multi_pos_weight_pair,ctrl_multi_pos_weight_pair,pos_treat_ctrl_multi_pair);
                 //info("after pair_treat_ctrl2");
                
                //output bedgraph

                write_multi_read_bedGraph_for_a_chromosome(peakIdx->get_chrom(i),pos_treat_ctrl_multi_pair, options.data_outdir, project_name);
                // info("after write bedgraph");
                
                double summit_pileup = 0.0;
                double summit_index = -1;
                double summit_pscore = 1.0;
                //info("pos_treat_ctrl_pair.pos size = " + std::to_string(pos_treat_ctrl_pair.pos.size()));
                
                if(pos_treat_ctrl_pair.pos.size() == 0 )
                {
                    continue;
                }
                
                int peak_start =pos_treat_ctrl_pair.pos[0];
                int peak_end = pos_treat_ctrl_pair.pos[pos_treat_ctrl_pair.pos.size()-1];
                
                for (size_t j=0; j< pos_treat_ctrl_pair.pos.size() ; j++)
                {
                    if(peak_start > pos_treat_ctrl_pair.pos[j])
                    {
                        peak_start = pos_treat_ctrl_pair.pos[j];
                    }
                    if(peak_end < pos_treat_ctrl_pair.pos[j])
                    {
                        peak_end = pos_treat_ctrl_pair.pos[j];
                    }
                    //get p-score
                    if(pos_treat_ctrl_pair.treat_v[j] > summit_pileup)
                    {
                        summit_index = j;
                        
                        summit_pileup = pos_treat_ctrl_pair.treat_v[j] * options.treat_scale;
                    }
                    
                }
                
                //info("after find summit ");
                double ctrl_expect = pos_treat_ctrl_pair.ctrl_v[summit_index] * options.ctrl_scale;
                //info("summit_pileup = " + std::to_string(int(summit_pileup)));
                //info("summit ctrl pileup = " + std::to_string(pos_treat_ctrl_pair.ctrl_v[summit_index]));
                //info("ctrl_expect = " + std::to_string(ctrl_expect));
                
                if(ctrl_expect < 0)
                {
                    info(" ctrl_expect less than 0");
                    
                    ctrl_expect = 0.0;
                }
                
                summit_pscore = -1.0 * log10_poisson_cdf(int(summit_pileup+0.5),roundTop((ctrl_expect+0.5),5),0);
                //info("summit_pscore = " + std::to_string(summit_pscore));
                
                if (p_chrom == peakIdx->get_chrom(i)){//check overlap
                
                    if (peak_start < p_end && peak_start > p_start ) //overlap with previous peak
                    {
                        p_start = std::min(p_start,peak_start);
                        p_end = std::max(p_end,peak_end);
                        p_ip_tag = p_ip_tag + (ip_tag + 0.1)* options.treat_scale;
                        p_input_tag = p_input_tag + (0.1 + input_tag)* options.ctrl_scale;
                        p_summit_pileup = p_summit_pileup + summit_pileup;
                        p_fe = std::max(fe,p_fe);
                        p_summit_pscore = std::max(summit_pscore,p_summit_pscore);
                        
                    }
                    else {
                        if(p_fe == 0.0) {p_fe = fe;}
                        ROUT << p_chrom << "\t" << p_start << "\t" << p_end << "\t" << p_ip_tag << "\t" << p_input_tag << "\t" << p_summit_pileup <<"\t" << p_fe << "\t" << p_summit_pscore << "\n";
                        
                        p_chrom = peakIdx->get_chrom(i);
                        p_start = peak_start;
                        p_end = peak_end;
                        p_ip_tag = (ip_tag + 0.1)* options.treat_scale;
                        p_input_tag = (0.1 + input_tag)* options.ctrl_scale;
                        p_summit_pileup = summit_pileup;
                        p_summit_pscore = summit_pscore;
                    }
                }
                else {
                    if(p_chrom != "") //not first peak
                    {
                        if(p_fe == 0.0) {p_fe = fe;}
                        
                        ROUT << p_chrom << "\t" << p_start << "\t" << p_end << "\t" << p_ip_tag << "\t" << p_input_tag << "\t" << p_summit_pileup <<"\t" << p_fe << "\t" << p_summit_pscore << "\n";
                    }
                        p_chrom = peakIdx->get_chrom(i);
                        p_start = peak_start;
                        p_end = peak_end;
                        p_ip_tag = (ip_tag + 0.1)* options.treat_scale;
                        p_input_tag = (0.1 + input_tag)* options.ctrl_scale;
                        p_summit_pileup = summit_pileup;
                        p_summit_pscore = summit_pscore;
                    
                }//end else overlap
                
            }//end pileup cutoff
            
            if(p_chrom != ""){
                if(p_fe == 0.0) {p_fe = fe;}
                ROUT << p_chrom << "\t" << p_start << "\t" << p_end << "\t" << p_ip_tag << "\t" << p_input_tag << "\t" << p_summit_pileup <<"\t" << p_fe << "\t" << p_summit_pscore << "\n";
                
            }
            p_chrom = "";
            p_start = 0;
            p_end = 0;
            p_summit_pscore = 0.0;
            p_ip_tag = 0;
            p_input_tag = 0;
            p_summit_pileup = 0;
            p_fe = 0.0;
            
        }//end all peaks
        
       // bedGraph_treat_f.close();
       // bedGraph_ctrl_f.close();
        ROUT.close();
        //f.close();
    }catch(std::ofstream::failure e ){
        std::cout << "Error in writing filtered peaks.\n" << std::endl;
    }

}
 

int run_EM_TEpeaks(opt_t options, ShortRead * treat, ShortRead * control,std::string peak_fname)
{

    std::vector<double> IP_peakReads;
    //std::vector<double> * IP_peakReads_ptr;
    std::vector<double> Input_peakReads;
    //std::vector<double> * Input_peakReads_ptr;
    std::vector<double> peak_reads_Prime;
    std::vector<double> ctrl_peak_reads_Prime;
    std::vector<int> peaks_with_multiReads;
    
    //mapping multi-reads to each peak
    MULTI_READ_MAP *  multi_read_mapTo_pids = new MULTI_READ_MAP();
    //unique reads start position
    std::map<int, std::vector<int> > * treat_peak_uniqReads_startPos = new std::map<int, std::vector<int> > ();
    std::map<int, std::vector<int> > * ctrl_peak_uniqReads_startPos = new std::map<int, std::vector<int> > ();
    
    std::map<int, peak_align_t> peak_multiAlign_map;
    
    info("Read in candidate peaks ... " + peak_fname);
    Candidate_Peaks * peakIdx = new Candidate_Peaks(peak_fname);

    
    info("Done reding candidate peaks.");
    //debug(" in run_EM_TEpeaks numofpeaks " + std::to_string(peakIdx->get_numofpeaks()));
    
    if (peakIdx == NULL) {
        std::cout << "error in read candidate peaks!" << std::endl;
        std::exit(1);
    }
    if (peakIdx->get_numofpeaks() < 10) {
        std::cout << "less than 10 candidate peak regions! Stop EM algorithm" << std::endl;
        std::exit(1);
    }
    //Candidate_Peaks peakIdx(peakbed);
    
    
    for(int i =0; i< peakIdx->get_numofpeaks() ; i++ )
    {
        IP_peakReads.push_back(0.0);
        Input_peakReads.push_back(0.0);
        peak_reads_Prime.push_back(0.0);
        ctrl_peak_reads_Prime.push_back(0.0);
        peaks_with_multiReads.push_back(0);
    }
    info("done peak read arrays...");
    //parse IP /input files for multi-reads
    //IP_peakReads_ptr = &IP_peakReads;
    
    
    read_distribution(treat, options.tfile,peakIdx,options,&IP_peakReads, multi_read_mapTo_pids,treat_peak_uniqReads_startPos);
    
    if (multi_read_mapTo_pids->size()==0) {
        
        delete peakIdx;
        
        delete multi_read_mapTo_pids;
        

        info("no multi_reads mapped to candidate peak regions.");
        return 0;
    }
    //EM re-distribute multi-reads
    std::vector<double> treat_multiAlgn_weight;
    std::vector<double> multiAlgn_final_weight;
    std::vector<int> multiAlgn_To_multiRead;
    std::vector<int> multiAlgn_To_peakID;
    std::vector<int> treat_multiAlgn_start; //keep start position
    
    //debug(" multi_read_mapTo_pids size " + std::to_string(multi_read_mapTo_pids->size()));
    
    prepare_EM_structure(peakIdx, IP_peakReads, multi_read_mapTo_pids, multiAlgn_To_multiRead, treat_multiAlgn_weight, multiAlgn_To_peakID,treat_multiAlgn_start);
    
    //save treat sample multiAlignments for each peak
    for(size_t k = 0; k <  multiAlgn_To_peakID.size(); k++)
    {
        int pid = multiAlgn_To_peakID[k];
        //std::cout << "#2 align_id = " << k << " pid = " << pid   << std::endl;
        
        if(peak_multiAlign_map.find(pid) != peak_multiAlign_map.end())
        {
            peak_multiAlign_map[pid].treat_multi_align_list.push_back(k);
        }
        else {
            peak_align_t pa;
            pa.treat_multi_align_list.push_back(k);
            peak_multiAlign_map.insert(std::pair<int,peak_align_t>(pid,pa));
        }
    }
    
    debug("after prepare EM structure numItr = " + std::to_string(options.numItr) );
    
    if (options.numItr == 0 ) { // eqauly distribute multi-alignment
        
        treat_multiAlgn_weight.clear();
        
        EqualWeight_reads(multiAlgn_To_multiRead, IP_peakReads, peak_reads_Prime, multiAlgn_To_peakID,treat_multiAlgn_weight);
    }
    else {
    EMestimate_read(multiAlgn_To_multiRead,multiAlgn_To_peakID,treat_multiAlgn_weight,IP_peakReads,peak_reads_Prime,options.numItr,peakIdx->peak_length);
        
        debug("after EMestimate_read");
    }

    int multiAlgn_first_idx = 0;
    
    for (size_t i=0; i< multiAlgn_To_multiRead.size(); i++) {
        //loop through all multi-alignments of each multi-read
        
        for (int j= multiAlgn_first_idx; j < multiAlgn_To_multiRead[i]; j++) {
            int pid = multiAlgn_To_peakID[j];
            
            
            peaks_with_multiReads[pid] = 1;
        }
        multiAlgn_first_idx = multiAlgn_To_multiRead[i];
    }
    
    info("Reding Input alignment file ...") ;
    //parse input file
    multi_read_mapTo_pids->clear();
    std::vector<double> ctrl_multiAlgn_weight;
    multiAlgn_final_weight.clear();
    multiAlgn_To_multiRead.clear();
    std::vector<int> ctrl_multiAlgn_start;
    multiAlgn_To_peakID.clear();
    //peak_reads_Prime.clear();
    
    //read in control sample
    read_distribution(control, options.cfile,peakIdx,options,&Input_peakReads,multi_read_mapTo_pids,ctrl_peak_uniqReads_startPos);
    
    debug("after control read distribution");
    
    //prepare EM structure for control sample
    prepare_EM_structure(peakIdx, Input_peakReads, multi_read_mapTo_pids, multiAlgn_To_multiRead, ctrl_multiAlgn_weight, multiAlgn_To_peakID,ctrl_multiAlgn_start);
    
    //info("before ctrl sample equal weighting");
    //do equal weighting only on control sample
    ctrl_multiAlgn_weight.clear();
    EqualWeight_reads(multiAlgn_To_multiRead, Input_peakReads, ctrl_peak_reads_Prime, multiAlgn_To_peakID,ctrl_multiAlgn_weight);
    
    //info("after ctrl sample equal weighting");
    
    //save control sample multiAlignments for each peak
    //info("before filter peaks");
    for(size_t k = 0; k <  multiAlgn_To_peakID.size(); k++)
    {
        int pid = multiAlgn_To_peakID[k];
        
        if(peak_multiAlign_map.find(pid) != peak_multiAlign_map.end())
        {
            peak_multiAlign_map[pid].ctrl_multi_align_list.push_back(k);
        }
        else {
            peak_align_t pa;
            pa.ctrl_multi_align_list.push_back(k);
            peak_multiAlign_map.insert(std::pair<int,peak_align_t>(pid,pa));
            
        }
    }
   // info("before filter peaks");
    //filter peaks
 filter_peaks(peak_multiAlign_map,treat_peak_uniqReads_startPos,ctrl_peak_uniqReads_startPos,peaks_with_multiReads,peakIdx, peak_reads_Prime, ctrl_peak_reads_Prime, treat_multiAlgn_weight, treat_multiAlgn_start, ctrl_multiAlgn_start,ctrl_multiAlgn_weight,options);
    
    
    //wirte bedgraph
    
    info("after filter peaks");
    delete peakIdx;
    delete multi_read_mapTo_pids;
    delete treat_peak_uniqReads_startPos;
    delete ctrl_peak_uniqReads_startPos;

    //write_filtered_peaks(peakIdx,data_outdir,project_name);
   // std::cout << "INFO  @ " << cur_time << "\tDone." << std::endl;
    
    return 1;
}

