//
//  CallerFromAlignments.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 5/31/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

//ifndef __TEToolkit_c____CallerFromAlignments__
//define __TEToolkit_c____CallerFromAlignments__

//include <stdio.h>
#include "ShortRead.h"
#include "PeakIO.h"

#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "Float64HashTable.h"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/unordered_map.hpp>
#include <thread>

//std::mutex peakMutex;



struct obs_exp_pair {
    obs_exp_pair(int observed,int expected) : tuple(observed,expected) {}
    
    boost::tuples::tuple<int,int> tuple;
};

struct pos_tc_t{
    std::vector<int> pos;
    std::vector<double> treat_v;
    std::vector<double> ctrl_v;
};

typedef boost::unordered_map< obs_exp_pair, double> ObsExpMap;

struct peak_content_t{
    int s;
    int e;
    double t_val;
    double c_val;
    int index;
};

struct call_peak_options_t
{
    std::vector<std::string> scoring_function_symbols;
    std::vector<double> score_cutoff_s;
    std::vector<double> lvl1_cutoff_s;
    std::vector<double> lvl2_cutoff_s;
    int lvl1_max_gap;
    int lvl2_max_gap;
    bool uniqOnly;
    int min_length;
    int max_gap;
    bool call_summits;
    bool auto_cutoff;
    int threads_to_use;

};

extern "C" {
    double roundTop(double num, int precision);
}

class CallerFromAlignments {

public:
    ShortRead *treat;            // FWTrack or PETrackI object for ChIP
    ShortRead *ctrl;             // FWTrack or PETrackI object for Control

    int  d ;                          // extension size for ChIP
    std::vector<int> ctrl_d_s ;                   // extension sizes for Control. Can be multiple values
    double treat_scaling_factor;       // scaling factor for ChIP
    std::vector<double> ctrl_scaling_factor_s;       // scaling factor for Control, corresponding to each extension size.
    double lambda_bg;                  // minimum local bias to fill missing values
    std::vector<std::string> chromosomes;                 // name of common chromosomes in ChIP and Control data
    std::vector<std::string> canoChromlist;
    double pseudocount;                // the pseudocount used to calcuate logLR, FE or logFE
    //std::string bedGraph_filename_prefix;     // prefix will be added to _pileup.bdg for treatment and _lambda.bdg for control

        //SHIFTCONTROL is obsolete
    int  end_shift ;                  // shift of cutting ends before extension
    bool trackline;                   // whether trackline should be saved in bedGraph
    bool save_bedGraph;               // whether to save pileup and local bias in bedGraph files
    bool save_SPMR ;                  // whether to save pileup normalized by sequencing depth in million reads
    bool no_lambda_flag;              // whether ignore local bias, and to use global bias instead
    bool PE_mode  ;                   // whether it's in PE mode, will be detected during initiation
        // temporary data buffer
    std::string chrom;                        // name of current chromosome
    pos_tc_t chr_pos_treat_ctrl;          // temporary [position, treat_pileup, ctrl_pileup] for a given chromosome
   // std::string bedGraph_treat_filename;
   // std::string bedGraph_control_filename;
    std::string  bedGraph_filename_prefix;
   // std::ofstream  bedGraph_treat_f;
   // std::ofstream  bedGraph_ctrl_f;
        //object bedGraph_treat            // file handler to write ChIP pileup
        //object bedGraph_ctrl             // file handler to write Control pileup
        // data needed to be pre-computed before peak calling
    Float64HashTable * pqtable ;                  // remember pvalue->qvalue convertion
    
    bool pvalue_all_done;             // whether the pvalue of whole genome is all calculated. If yes, it's OK to calculate q-value.

    std::map<double,int> pvalue_npeaks;               // record for each pvalue cutoff, how many peaks can be called
    std::map<double,int> pvalue_length;               // record for each pvalue cutoff, the total length of called peaks
    
    double optimal_p_cutoff;           // automatically decide the p-value cutoff ( can be translated into qvalue cutoff ) based
                                         // on p-value to total peak length analysis. 
   // std::string cutoff_analysis_filename;     // file to save the pvalue-npeaks-totallength table

    double test_time;
    
    std::map<std::string,std::string> pileup_data_files;           // Record the names of temporary files for storing pileup values of each chromosome
    
    
    CallerFromAlignments(std::vector<std::string> genome_chromlist, ShortRead * treat,
                         ShortRead * ctrl,
                         std::vector<int> ctrl_d_s,
                          std::vector<double> ctrl_scaling_factor_s,
                         int d = 200,
                         double treat_scaling_factor = 1.0,
                         double pseudocount = 1.0,
                         int end_shift = 0,
                         double lambda_bg = 0,
                         bool no_lambda_flag = false,
                         bool save_bedGraph = true,
                         std::string  bedGraph_filename_prefix = "PREFIX");
    
    ~CallerFromAlignments();
    
    void call_peaks(PeakIO * peaks, call_peak_options_t copts);
    void call_broadpeaks ( PeakIO *peaks, call_peak_options_t copts );
private:
    ObsExpMap pscore_map;
    
    void merge_peaklist(PeakIO *peaks, std::vector<PeakIO *> peaks_list,bool uniqOnly);

    double get_pscore ( int observed, double expectation );
    double get_pscore_v2 ( int observed, double expectation );
    std::vector<double> __cal_qscore ( std::vector<double> array1, std::vector<double> array2 );
    
    void __pre_computes ( int max_gap = 50, int min_length = 200 );
    void __cal_pvalue_qvalue_table();
    std::vector<double> __cal_FE ( std::vector<double> array1, std::vector<double> array2);
    std::vector<double> __cal_pscore ( std::vector<double> array1, std::vector<double> array2, bool uniqOnly);
    
    void __pileup_treat_ctrl_a_chromosome(std::string chrom,bool uniqOnly);
    pos_tc_t __pileup_treat_ctrl_a_chromosome_multiThread(std::string chrom, bool uniqOnly);
    
    void __chrom_pair_treat_ctrl_multiThread(std::pair<std::vector<int> , std::vector<double> > treat_pv, std::pair<std::vector<int>, std::vector<double> > ctrl_pv, pos_tc_t &res_pos_tc);
    
    void __chromlist_call_peak_using_certain_criteria(std::vector<std::string> chromlist,PeakIO *peaks, call_peak_options_t copts );
    
    void __chrom_pair_treat_ctrl ( std::pair<std::vector<int>,std::vector<double> > treat_pv, std::pair<std::vector<int>,std::vector<double> > ctrl_pv );
    
    void __chrom_call_peak_using_certain_criteria (PeakIO  *peaks, std::string chrom, call_peak_options_t copts, bool save_bedGraph );

    
    bool __close_peak_with_subpeaks (std::vector<peak_content_t> peak_content,PeakIO *peaks, int min_length, std::string chrom, int smoothlen, std::vector<std::vector<double> > score_array_s, std::vector<double> score_cutoff_s,float min_valley = 0.9 );
    
    bool __close_peak_wo_subpeaks (bool uniqOnly,std::vector<peak_content_t> peak_content, PeakIO *peaks, int min_length, std::string chrom, int smoothlen, std::vector<std::vector<double> > score_array_s, std::vector<double> score_cutoff_s);
    
    
    void __chrom_call_broadpeak_using_certain_criteria (PeakIO  *peaks,PeakIO  *peaks2, std::string chrom, call_peak_options_t copts, bool save_bedGraph );
    void __chromlist_call_broadpeak_using_certain_criteria(std::vector<std::string> chromlist,PeakIO *peaks, call_peak_options_t copts );
    bool __close_peak_for_broad_region (bool uniqOnly,std::vector<peak_content_t> peak_content, PeakIO *peaks, int min_length, std::string chrom, int smoothlen, std::vector<std::vector<double> > score_array_s);
    double mean_from_value_length(std::vector<double> array1, std::vector<int> lenlist);
    void __write_bedGraph_for_a_chromosome(std::string chrom,pos_tc_t chr_pos_treat_ctrl_multiThread);
    void __add_broadpeak (PeakIO *peaks, std::string chrom, PeakContent lvl2peak, std::vector<PeakContent> lvl1peakset);
};
//endif /* defined(__TEToolkit_c____CallerFromAlignments__) */
