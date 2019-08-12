//
//  PeakDetect.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 5/3/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//
// This code is modified from MACS2.

#include <math.h>
#include "PeakDetect.h"
#include "myLog.h"
#include "PeakModel.h"
#include "assert.h"
#include "CallerFromAlignments.h"
#include <sstream>
#include <iostream>
#include <iomanip>

    //Class to do the peak calling.

PeakDetect::PeakDetect(ShortRead *treat, ShortRead * control, opt_t *opt)
{
    this->opt = opt;

    this->treat = treat;
    this->control = control;
    this->ratio_treat2control = 0.0;
    this->peaks = new PeakIO();
    
    //this->PE_MODE = opt->PE_mode;
    //this->scoretrack = NULL;
    
    //this.log_pvalue = opt->opt->log_pvalue;    // -log10pvalue
   // this.opt->log_qvalue = opt->opt->log_qvalue;    // -log10qvalue
    
    this->d = opt->d;
    //this.gsize = opt->gsize;
    //this.end_shift = opt->shift;
    //this.nolambda = opt->nolambda;
    //this.opt->smalllocal = opt->smalllocal;
    //this.opt->largelocal = opt->largelocal;
    
    //ADD if (opt->nolambda){
    //ADD    info("#3 !!!! DYNAMIC LAMBDA IS DISABLED !!!!");
    //ADD }
}

PeakDetect::~PeakDetect(){
    delete peaks;
}

//find candidate peak regions for considering multi-reads
void PeakDetect::candidate_peakregions()
{

    if (this->peaks != nullptr) {
        delete this->peaks;
        this->peaks = new PeakIO();
    }
    
    __call_peaks_w_control(false);

}

void PeakDetect::call_peaks ()
{
        /*Call peaks function.

        Scan the whole genome for peaks. RESULTS WILL BE SAVED IN
        final_peaks and final_negative_peaks.
        */
    
//ADD        if(control)                // w/ control
//ADD        {
//ADD            __call_peaks_w_control ();
 //ADD       }
 //ADD       else{                           // w/o control
 //ADD           __call_peaks_wo_control ();
 //ADD       }
    __call_peaks_w_control (true);
}

void PeakDetect::__call_peaks_w_control (bool uniqOnly)
{
        /*"""To call peaks with control data.
        
        While calculating pvalue:

        First, t and c will be adjusted by the ratio between total
        reads in treatment and total reads in control, depending on
        --to-small option.

        Then, t and c will be multiplied by the smallest peak size --
        d.

        Finally, a poisson CDF is applied to calculate one-side pvalue
        for enrichment.
        """*/
    
    //int i;
    double lambda_bg; //, effective_depth_in_million;
    double treat_scale, d;
    std::vector<double> ctrl_scale_s;
    std::vector<int> ctrl_d_s;
    double  control_total;
    double treat_sum;              // approx sum of treatment pileup values
    double control_sum;            // approx sum of control pileup values
    bool save_bedGraph;
    
    
            
    //treat_total   = this->treat->get_total();
    
    if (opt->PE_mode){
        
        d = treat->fraglength;
        control_total = this->control->get_total() * 2; // in PE mode, entire fragment is counted as 1
                                                   // in treatment whereas both ends of fragment are counted in control/input. ??
        treat_sum = treat->length;
        control_sum = control_total * treat->average_template_length;
        ratio_treat2control = double(treat_sum)/control_sum;
    }
    else { //SE
        d = opt->d;
        control_total = control->get_total();
        treat_sum = treat->get_total() * d;
        control_sum = control->get_total() * d;
        ratio_treat2control = double(treat_sum)/control_sum;
    }

    if (opt->ratio != 1.0){
        ratio_treat2control = opt->ratio;
    }

    if (opt->tocontrol){
            // if MACS decides to scale treatment to control data because treatment is bigger
           // effective_depth_in_million = control_total / 1000000.0;
        lambda_bg = double( control_sum )/ opt->gsize;
            treat_scale = 1.0/ratio_treat2control;
        
        opt->treat_scale = treat_scale;
        opt->ctrl_scale = 1.0;
    }
    else{
            // if MACS decides to scale control to treatment because control sample is bigger
       // effective_depth_in_million = treat_total / 1000000.0;
        lambda_bg = double( treat_sum )/ opt->gsize;
        treat_scale = 1.0;
        
        opt->treat_scale = 1.0;
        opt->ctrl_scale = ratio_treat2control;
    }
    // prepare d_s for control data
    if (opt->smalllocal != 0 ){
        assert(d <= opt->smalllocal); //, "slocal can't be smaller than d!"
    }
    if (opt->largelocal != 0){
        assert(d <= opt->largelocal); // , "llocal can't be smaller than d!"
        assert(opt->smalllocal <= opt->largelocal); //, "llocal can't be smaller than slocal!"
    }

    // Now prepare a list of extension sizes
    //std::vector<int> ctrl_d_s; // = [ d ];   // note, d doesn't make sense in PE mode.
    ctrl_d_s.push_back(d);
    // And a list of scaling factors for control
    //std::vector<double> ctrl_scale_s;
        // d
    double tmp_v = 1.0;
    if (! opt->tocontrol){
            // if user wants to scale everything to ChIP data
        tmp_v = ratio_treat2control;
    }

    ctrl_scale_s.push_back( tmp_v );

        // slocal size local
    if (opt->smalllocal != 0){
        ctrl_d_s.push_back( opt->smalllocal );
        
        if (! opt->tocontrol){
                // if user want to scale everything to ChIP data
            tmp_v = 1.0 * d/opt->smalllocal * ratio_treat2control;
        }
        else{
            tmp_v = 1.0 * d/opt->smalllocal;
            
        }
        ctrl_scale_s.push_back( tmp_v );
    }

        // llocal size local
    if (opt->largelocal != 0 and opt->largelocal > opt->smalllocal){
        ctrl_d_s.push_back( opt->largelocal );
        if (not opt->tocontrol){
                // if user want to scale everything to ChIP data
            tmp_v = 1.0 * d /opt->largelocal * ratio_treat2control;
        }
        else{
            tmp_v = 1.0 * d/opt->largelocal;
        }
        ctrl_scale_s.push_back( tmp_v );
    }
    
//ADD    if (opt->nolambda){
//ADD        ctrl_d_s.clear();// = [];
//ADD        ctrl_scale_s.clear();
//ADD    }
    info("opt.gsize = " + std::to_string(opt->gsize));

    info("lambda_bg " + std::to_string(lambda_bg) + " treat_sum = " + std::to_string(treat_sum));
    std::string bdg_prefix = opt->data_outdir + "/" + opt->project_name;
    if(uniqOnly)
    {
        save_bedGraph = true;
        
    }
    else {
        save_bedGraph = false;
    }
        
      CallerFromAlignments  scorecalculator( opt->chromlist, treat, control,ctrl_d_s,ctrl_scale_s,
                                                d,
                                                treat_scale,
                                            
                                                1.0,
                                                opt->shift,
                                            lambda_bg, false, save_bedGraph, bdg_prefix);


       // if opt->trackline: scorecalculator.enable_trackline()
        // call peaks
//ADD    if (opt->call_summits){
//ADD        info("#3 Going to call summits inside each peak ...");
//ADD    }
    std::vector<std::string> scoring_function_symbols;
    std::vector<double> cutoff_list;
    
    call_peak_options_t call_peak_opts;
    
    call_peak_opts.min_length = d;
    call_peak_opts.max_gap = opt->tsize;
    call_peak_opts.call_summits = false;
    call_peak_opts.auto_cutoff = false;
    call_peak_opts.uniqOnly = uniqOnly;
    call_peak_opts.threads_to_use = opt->threadNum;
    
    if (opt->log_qvalue && uniqOnly ){
        
        info("#4 Call peaks with given  q-value cutoff : " + std::to_string (opt->fdr));
      
        scoring_function_symbols.push_back("q");
        cutoff_list.push_back(-1.0*log10(opt->fdr));
        
    }
    else{
        
        
        stringstream stream;
        
        stream << std::fixed << std::setprecision(10) << opt->pval;
        
        string s = stream.str();
        info("#4 Call peaks with given p-value cutoff: " + s);
        
        scoring_function_symbols.push_back("p");
        cutoff_list.push_back(-1.0 * log10(opt->pval));

    }

    call_peak_opts.scoring_function_symbols = scoring_function_symbols;
    call_peak_opts.score_cutoff_s = cutoff_list;
    

    if(opt->broad)
    {
        call_peak_opts.lvl1_cutoff_s = cutoff_list;
        call_peak_opts.lvl1_max_gap=opt->tsize;
        call_peak_opts.lvl2_max_gap=d * 4;
        call_peak_opts.lvl2_cutoff_s.push_back(opt->log_broadcutoff);
        
        scorecalculator.call_broadpeaks(this->peaks,call_peak_opts);
    }
    else {
        scorecalculator.call_peaks(this->peaks,call_peak_opts);
    }
    
    //scorecalculator.destroy();
    
}
/*
void PeakDetect::__call_peaks_wo_control (bool uniqOnly)
{
    /* """To call peaks without control data.
     A peak info type is a: dictionary
     
     key value: chromosome
     
     items: (peak start,peak end, peak length, peak summit, peak
     height, number of tags in peak region, peak pvalue, peak
     fold_enrichment) <-- tuple type
     
     While calculating pvalue:
     
     First, t and c will be adjusted by the ratio between total
     reads in treatment and total reads in control, depending on
     --to-small option.
     
     Then, t and c will be multiplied by the smallest peak size --
     d.
     
     Finally, a poisson CDF is applied to calculate one-side pvalue
     for enrichment.
     """ */
/*
    double lambda_bg, effective_depth_in_million;
    
    double d;
    
    std::vector<double> ctrl_scale_s, ctrl_d_s;
    
    if (this.PE_MODE){
        d = 0;
    }
    else{
        d = this.d;
    }
    int treat_length = treat->length;
    int treat_total = treat->get_total();
    
    effective_depth_in_million = treat_total / 1000000.0;
    
    // global lambda
    if (this.PE_MODE){
        //    # this an estimator, we should maybe test it for accuracy?
        lambda_bg = treat_length / opt->gsize;
    }
    else{
        lambda_bg = 1.0 * d * treat_total / opt->gsize;
    }
    treat_scale = 1.0;
    
    // slocal and d-size local bias are not calculated!
    // nothing done here. should this match w control??
    
    if (not opt->nolambda){
        if (PE_MODE){
            ctrl_scale_s.push_back(1.0 * treat_length/ (opt->largelocal * treat_total * 2));
        }
        else {
            ctrl_scale_s.push_back(1.0 * d / opt->largelocal);
        }
        ctrl_d_s.push_back(opt->largelocal);
    }
    //calculate peak scores, do each chromosome separately (can be parallellized)
    CallerFromAlignments scorecalculator ( treat, NULL,
                                          d = d, ctrl_d_s = ctrl_d_s,
                                          treat_scaling_factor = 1.0,
                                          ctrl_scaling_factor_s = ctrl_scale_s,
                                          pseudocount = 1.0,
                                          end_shift = opt->shift,
                                          lambda_bg = lambda_bg,
                                          no_lambda_flag = false,
                                          save_bedGraph = opt->store_bdg,
                                          bedGraph_filename_prefix = opt->name,
                                          bedGraph_treat_filename = opt->bdg_treat,
                                          bedGraph_control_filename = opt->bdg_control,
                                          cutoff_analysis_filename = opt->cutoff_analysis_file,
                                          save_SPMR = opt->do_SPMR
                                          );
    
    //if opt->trackline: scorecalculator.enable_trackline()
    // call peaks
    if(opt->call_summits){
        info("#3 Going to call summits inside each peak ...");
    }
    
    std::vector<std::string> scoring_function_symbols;
    std::vector<double> cutoff_list;
    
    if (opt->log_pvalue){
        /* if (opt->broad){
         info("#3 Call broad peaks with given level1 -log10pvalue cutoff and level2: %.5f, %.5f..." % (opt->log_pvalue,opt->log_broadcutoff) );
         peaks = scorecalculator.call_broadpeaks(['p',], lvl1_cutoff_s=[opt->log_pvalue,],lvl2_cutoff_s=[opt->log_broadcutoff,],min_length=d,
         lvl1_max_gap=opt->tsize,lvl2_max_gap=d*4);
         }*/
        //else{
/*        info("#3 Call peaks with given -log10pvalue cutoff: %.5f ..." % opt->log_pvalue);
        scoring_function_symbols.push_back("p");
        cutoff_list.push_back(opt->log_pvalue);
        //}
    }
    
    
    if (opt->log_qvalue){
        /*if (opt->broad){
         info("#3 Call broad peaks with given level1 -log10qvalue cutoff and level2: %f, %f..." % (opt->log_qvalue,opt->log_broadcutoff) );
         peaks = scorecalculator.call_broadpeaks(['q',], lvl1_cutoff_s=[opt->log_qvalue,],lvl2_cutoff_s=[opt->log_broadcutoff,],min_length=d,
         lvl1_max_gap=opt->tsize,lvl2_max_gap=d*4);
         }
         else{*/
        
        //}
/*        scoring_function_symbols.push_back("q");
        cutoff_list.push_back(opt->log_qvalue);
        
    }
    scorecalculator.call_peaks(peaks, scoring_function_symbols, cutoff_list,
                               min_length=d,
                               max_gap=opt->tsize,
                               call_summits=opt->call_summits,
                               auto_cutoff=opt->cutoff_analysis );
    
    
}
*/
