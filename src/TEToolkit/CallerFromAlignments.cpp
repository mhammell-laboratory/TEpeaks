//
//  CallerFromAlignments.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 5/31/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

//Modified from MACS2

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include <fstream>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <map>
#include <vector>
#include <functional>
#include <boost/lexical_cast.hpp>
#include "CallerFromAlignments.h"
#include "myLog.h"
#include "zeroin.h"
#include "cStatistics.h"
#include "Constants.h"
#include <sys/time.h>
#include <thread>



using boost::lexical_cast;

//typedef std::numeric_limits< double > dbl;


double roundTop(double num, int precision)
{
    
    return floor(num * pow(10.0f,precision) + .5f)/pow(10.0f,precision);
}

bool operator==(const obs_exp_pair &a , const obs_exp_pair &b)
{
    
    return a.tuple.get<0>() == b.tuple.get<0>() &&
    a.tuple.get<1>() == b.tuple.get<1>();
    
}

std::size_t hash_value(const obs_exp_pair &a)
{
    std::size_t seed =0;
    boost::hash_combine(seed,a.tuple.get<0>());
    boost::hash_combine(seed,a.tuple.get<1>());
    
    return seed;
}



void apply_multiple_cutoffs ( std::vector<std::vector<double>> multiple_score_arrays, std::vector<double> multiple_cutoffs,std::vector<int> &indices)
{
    int array_length = multiple_score_arrays[0].size();
    std::vector<int> tmp;
    
    for (int i=0; i < array_length; i++) {
        tmp.push_back(0);
    }
    
    for (int i= 0; i < array_length;i++)
    {
        
        for (size_t j=0;j < multiple_cutoffs.size(); j++){
            if (multiple_score_arrays[j][i] > multiple_cutoffs[j]) {
                tmp[i] = 1;
                break;
            }
        }
        
    }
    
    for (int i =0 ; i < array_length; i++) {
        if (tmp[i] == 1) {
            indices.push_back(i);
        }
    }
}

std::vector<int> apply_FE_cutoff ( std::vector<double> array1, std::vector<double> array2 ,std::vector<int> indices)
{
    std::vector<int> ret;
    if (array1.size() != array2.size()) {
        error("Two arraies shouble have the same length.");
        std::exit(1);
    }
    
    for (size_t i=0 ; i< indices.size();i++)
    {
         double s = 1.0 * (array1[indices[i]]+1)/(array2[indices[i]]+1);
         if( s > 50)
         {
              ret.push_back(indices[i]);
         }
    }
    return ret;
}
double CallerFromAlignments::get_pscore_v2 ( int observed, double expectation )
{
    /*"""Get p-value score from Poisson test. First check existing
     table, if failed, call poisson_cdf function, then store the result
     in table.
     
     """*/
    
    double score;
    
    score = -1.0 * log10_poisson_cdf(observed,roundTop(expectation,5),0);
    
    return score;
}
double CallerFromAlignments::get_pscore ( int observed, double expectation )
{
    /*"""Get p-value score from Poisson test. First check existing
     table, if failed, call poisson_cdf function, then store the result
     in table.
     
     """*/
    
    double score;
    //long key_value;
   // double expectation_round = round(expectation * 100000)/100000;
   // std::string key_str = lexical_cast<std::string>(observed) + ":" + lexical_cast<std::string>(expectation_round);
    
    //std::hash<std::string> str_hash;
    //size_t key_value =  str_hash( key_str );
    //info("get_pscore");
      
    //obs_exp_pair key (observed,roundTop(expectation,5));
    
    int e = floor((expectation * pow(10.0f,5) + .5f)/10);
    
    //std::cout << std::setprecision(5) << expectation << "\t e = " << e << std::endl;
    
    obs_exp_pair key (observed,e);
    
    obs_exp_pair key1 (observed,e+1);
    obs_exp_pair key2 ( observed, e-1);
   
    if (pscore_map.find(key) != pscore_map.end()) {
        return pscore_map[key];
    }
    if (pscore_map.find(key1) != pscore_map.end()) {
        return pscore_map[key1];
    }
    if (pscore_map.find(key2) != pscore_map.end()) {
        return pscore_map[key2];
    }
    else {
        
        score = -1.0 * log10_poisson_cdf(observed,roundTop(expectation,5),0);
        pscore_map.emplace(key,score) ; 
    }
    
    return score;
}

CallerFromAlignments::CallerFromAlignments(std::vector<std::string> genome_chromlist,ShortRead * treat, ShortRead * ctrl,std::vector<int> ctrl_d_s,
                                           std::vector<double> ctrl_scaling_factor_s,
                                           int d ,
                                           double treat_scaling_factor ,
                                           double pseudocount ,
                                           int end_shift,
                                           double lambda_bg , bool no_lambda_flag,
                                           bool save_bedGraph,
                                           std::string  bedGraph_filename_prefix )
{
        /*"""Initialize.

        A calculator is unique to each comparison of treat and
        control. Treat_depth and ctrl_depth should not be changed
        during calculation.

        treat and ctrl are two ShortRead objects.

        treat_depth and ctrl_depth are effective depth in million:
                                    sequencing depth in million after
                                    duplicates being filtered. If
                                    treatment is scaled down to
                                    control sample size, then this
                                    should be control sample size in
                                    million. And vice versa.

        d, sregion, lregion: d is the fragment size, sregion is the
                             small region size, lregion is the large
                             region size
                                    
        pseudocount: a pseudocount used to calculate logLR, FE or
                     logFE. Please note this value will not be changed
                     with normalization method. So if you really want
                     to set pseudocount 1 per million reads, set it
                     after you normalize treat and control by million
                     reads by `change_normalizetion_method(ord('M'))`.

        """ */
    if (ctrl_d_s.size() == 0) {
        ctrl_d_s.push_back(200);
        ctrl_d_s.push_back(1000);
        ctrl_d_s.push_back(10000);
    }
    if (ctrl_scaling_factor_s.size() == 0) {
        ctrl_scaling_factor_s.push_back(1.0);
        ctrl_scaling_factor_s.push_back(0.2);
        ctrl_scaling_factor_s.push_back(0.02);
    }
    
    
    if (treat->isPE) {
        PE_mode = true;
    }
    else {
        PE_mode = false;
    }
    this->treat = treat;
    this->ctrl = ctrl;
    
    if (this->ctrl == NULL) {
        this->ctrl = this->treat;
    }
        
    this->trackline = false;
    this->d = d;
    this->ctrl_d_s = ctrl_d_s;

    this->treat_scaling_factor = treat_scaling_factor;
    this->ctrl_scaling_factor_s= ctrl_scaling_factor_s;
    this->end_shift = end_shift;
    this->lambda_bg = lambda_bg;
    this->pqtable = nullptr;
    
    this->save_bedGraph = save_bedGraph;
    this->save_SPMR = save_SPMR;
    this->bedGraph_filename_prefix =  bedGraph_filename_prefix;
    //tmp_bytes = bedGraph_treat_filename.encode('UTF-8')
    //print bedGraph_treat_filename, tmp_bytes
    //this->bedGraph_treat_filename = bedGraph_treat_filename;
    //tmp_bytes = bedGraph_control_filename.encode('UTF-8')
    //print bedGraph_control_filename, tmp_bytes
    //this->bedGraph_control_filename = bedGraph_control_filename;

    this->no_lambda_flag = no_lambda_flag;

    this->pseudocount = pseudocount;

    std::vector<std::string> chr1 = this->treat->get_chrom_names();
    std::vector<std::string> chr2 = this->ctrl->get_chrom_names();
    
    
    for (auto ch : chr1) {
        //if (chr2.find(ch) != chr2.end())
        if (std::find(chr2.begin(), chr2.end(), ch) != chr2.end() )
        {
            this->chromosomes.push_back(ch);
            if (std::find(genome_chromlist.begin(),genome_chromlist.end(),ch) != genome_chromlist.end()) {
                this->canoChromlist.push_back(ch);
            }
            
           // if(ch.length() < 6)
           // {
           //     this->canoChromlist.push_back(ch);
           // }
        }
    }

    this->test_time = 0;

    double f = 0.3;
    // step for optimal cutoff is 0.3 in -log10pvalue, we try from pvalue 1E-10 (-10logp=10) to 0.5 (-10logp=0.3)
    while(f < 10) {
        this->pvalue_length.insert(std::pair<double,int>(f,0));
        this->pvalue_npeaks.insert(std::pair<double,int>(f,0));
        f += 0.3;
    }
    this->optimal_p_cutoff = 0;
    //this->cutoff_analysis_filename = cutoff_analysis_filename;
}

CallerFromAlignments::~CallerFromAlignments(){
        /*"""Remove temparary files for pileup values of each chromosome.

        Note: This function MUST be called if the class object won't
        be used anymore.

        """*/
    
    //std::string f;

    for(auto f : pileup_data_files)
    {
        std::ifstream ff(f.second);
        if( ff.good() ){
            std::remove(f.second.c_str());
        }
    }
    
    delete pqtable;
    
    return ;
}

void CallerFromAlignments::__chrom_pair_treat_ctrl ( std::pair<std::vector<int>,std::vector<double> > treat_pv, std::pair<std::vector<int>,std::vector<double> > ctrl_pv )
{
        /*"""*private* Pair treat and ctrl pileup for each region.
        
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
            chr_pos_treat_ctrl.pos.push_back(treat_pv.first[it]);
            chr_pos_treat_ctrl.treat_v.push_back(treat_pv.second[it]);
            chr_pos_treat_ctrl.ctrl_v.push_back(ctrl_pv.second[ic]);

                // call for the next p1 and v1
            it += 1;
        }
        else {
        if(treat_pv.first[it] > ctrl_pv.first[ic]){
                // clip a region from pre_p to p2, then set pre_p as p2.
            
            chr_pos_treat_ctrl.pos.push_back(ctrl_pv.first[ic]);
            chr_pos_treat_ctrl.treat_v.push_back(treat_pv.second[it]);
            chr_pos_treat_ctrl.ctrl_v.push_back(ctrl_pv.second[ic]);

                // call for the next p2 and v2
                ic += 1;
        }
        else{
                // from pre_p to p1 or p2, then set pre_p as p1 or p2.
            chr_pos_treat_ctrl.pos.push_back(treat_pv.first[it]);
            chr_pos_treat_ctrl.treat_v.push_back(treat_pv.second[it]);
            chr_pos_treat_ctrl.ctrl_v.push_back(ctrl_pv.second[ic]);

                // call for the next p1, v1, p2, v2.
            it += 1;
            ic += 1;

        }
        }
    }

}

void CallerFromAlignments::__chrom_pair_treat_ctrl_multiThread ( std::pair<std::vector<int>,std::vector<double> > treat_pv, std::pair<std::vector<int>,std::vector<double> > ctrl_pv,pos_tc_t &res_pos_tc  )
{
    /*"""*private* Pair treat and ctrl pileup for each region.
     
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
            res_pos_tc.pos.push_back(treat_pv.first[it]);
            res_pos_tc.treat_v.push_back(treat_pv.second[it]);
            res_pos_tc.ctrl_v.push_back(ctrl_pv.second[ic]);
            
            // call for the next p1 and v1
            it += 1;
        }
        else {
            if(treat_pv.first[it] > ctrl_pv.first[ic]){
                // clip a region from pre_p to p2, then set pre_p as p2.
                
                res_pos_tc.pos.push_back(ctrl_pv.first[ic]);
                res_pos_tc.treat_v.push_back(treat_pv.second[it]);
                res_pos_tc.ctrl_v.push_back(ctrl_pv.second[ic]);
                
                // call for the next p2 and v2
                ic += 1;
            }
            else{
                // from pre_p to p1 or p2, then set pre_p as p1 or p2.
                res_pos_tc.pos.push_back(treat_pv.first[it]);
                res_pos_tc.treat_v.push_back(treat_pv.second[it]);
                res_pos_tc.ctrl_v.push_back(ctrl_pv.second[ic]);

                it += 1;
                ic += 1;
                
            }
        }
    }
    
}

void CallerFromAlignments::__pileup_treat_ctrl_a_chromosome ( std::string chrom , bool uniqOnly )
{
        /*"""After this function is called, chr_pos_treat_ctrl will
        be reset and assigned to the pileup values of the given
        chromosome.
        
        """*/
    
    std::pair<std::vector<int>, std::vector<double> > treat_pv, ctrl_pv;
    //long i;
    //float t;
    std::ifstream f;

    if (std::find(chromosomes.begin(),chromosomes.end(),chrom) == chromosomes.end()) {
        error("chromosome is not valid : " + chrom);
        std::exit(1);
    }
    //assert chrom in self.chromosomes, "chromosome %s is not valid." % chrom;

    // check backup file of pileup values. If not exists, create
    // it. Otherwise, load them instead of calculating new pileup values.
    if (uniqOnly) {
    
    if (pileup_data_files.find( chrom ) != pileup_data_files.end())
    {
        
     try{
         f.open( pileup_data_files[ chrom ],std::ifstream::in );
         chr_pos_treat_ctrl.pos.clear();
         chr_pos_treat_ctrl.treat_v.clear();
         chr_pos_treat_ctrl.ctrl_v.clear();
         
        while(! f.eof()){
            std::string line, ss_pos,ss_treat,ss_ctrl;
            std::stringstream ss;
            
            if (! std::getline(f,line)){
                break;
            }
            
            ss << line;
            std::getline(ss,ss_pos,'\t');
            std::getline(ss,ss_treat,'\t');
            std::getline(ss,ss_ctrl,'\t');
            
            chr_pos_treat_ctrl.pos.push_back(std::stol(ss_pos));
            chr_pos_treat_ctrl.treat_v.push_back(std::stod(ss_treat));
            chr_pos_treat_ctrl.ctrl_v.push_back(std::stod(ss_ctrl));
            
         }
            
         f.close();
            return;
        }
        catch(std::ifstream::failure e ){
            struct timeval start;
            gettimeofday(&start, NULL);
            std::string temp_filename = lexical_cast<std::string>(start.tv_sec) + "." + lexical_cast<std::string> ( start.tv_usec);
            //temp_fd, temp_filename = mkstemp();
            //os.close(temp_fd);
            pileup_data_files.insert(std::pair<std::string,std::string>(chrom,temp_filename));
        }
    }
    else{
        //temp_fd, temp_filename = mkstemp();
        //os.close(temp_fd);
        //pileup_data_files[ chrom ] = temp_filename;
        struct timeval start;
        gettimeofday(&start, NULL);
        std::string temp_filename = lexical_cast<std::string>(start.tv_sec) + "." + lexical_cast<std::string> ( start.tv_usec);
        pileup_data_files.insert(std::pair<std::string,std::string>(chrom,temp_filename));
       
    }
    }

    // reset or clean existing chr_pos_treat_ctrl
    chr_pos_treat_ctrl.pos.clear();
    chr_pos_treat_ctrl.treat_v.clear();
    chr_pos_treat_ctrl.ctrl_v.clear();

    std::vector<double> scale_factors;
    std::vector<int> ds;
    ds.push_back(this->d);
    scale_factors.push_back(treat_scaling_factor);
    
    if (PE_mode){
        
        treat_pv = treat->pileup_a_chromosome_pe ( chrom, scale_factors, 0.0 );
        debug("treat_pv size = " + std::to_string(treat_pv.first.size()));
    }
    else{
        treat_pv = treat->pileup_a_chromosome( chrom, ds, scale_factors,uniqOnly, 0.0,true, this->end_shift );
    }
    
    if (not no_lambda_flag){
        
            if (PE_mode){
                // note, we pileup up PE control as SE control because
                // we assume the bias only can be captured at the
                // surrounding regions of cutting sites from control experiments.
                
                ctrl_pv = ctrl->pileup_a_chromosome_c( chrom, this->ctrl_d_s, this->ctrl_scaling_factor_s, this->lambda_bg );
                debug("ctrl_pv size = " + std::to_string(ctrl_pv.first.size()));
                for (size_t i=0 ; i < ctrl_pv.first.size(); i ++ ) {
                    
                    std::cout  << ctrl_pv.first[i] << "\t pos " << std::endl;
                }
            }
            else{
                
                ctrl_pv = ctrl->pileup_a_chromosome( chrom, this->ctrl_d_s, this->ctrl_scaling_factor_s,uniqOnly,
                                                         this->lambda_bg,
                                                        false );
                
            }
    }
    else{
        std::vector<int> pos;
        pos.push_back(treat_pv.first[treat_pv.first.size()-1]);
        std::vector<double> val;
        val.push_back(this->lambda_bg);
        
        ctrl_pv = std::pair<std::vector<int>,std::vector<double> > (pos,val);
        
        
    }

    __chrom_pair_treat_ctrl( treat_pv, ctrl_pv);
    
    // clean treat_pv and ctrl_pv
    treat_pv.first.clear(); // = [];
    treat_pv.second.clear();
    ctrl_pv.first.clear(); //  = [];
    ctrl_pv.second.clear();
    
    info( "chr_pos_treat_ctrl pos size = " + std::to_string(chr_pos_treat_ctrl.pos.size()) );
    // save data to temporary file
    if (uniqOnly) {
    
    try{
        std::ofstream f(this->pileup_data_files[ chrom ],std::ofstream::out);
        //    cPickle.dump( chr_pos_treat_ctrl, f , protocol=2 );
        
        for (size_t i=0 ; i < chr_pos_treat_ctrl.pos.size(); i ++ ) {
            
            //double tr_v = round(chr_pos_treat_ctrl.treat_v[i] * 100000)/100000;
            //double ct_v = round(chr_pos_treat_ctrl.ctrl_v[i] * 100000)/100000;
            
            double tr_v = chr_pos_treat_ctrl.treat_v[i];
            double ct_v = chr_pos_treat_ctrl.ctrl_v[i];
                                
            f  << chr_pos_treat_ctrl.pos[i] << "\t" << tr_v << "\t" << std::setprecision(5) << ct_v << "\n";
        }
        f.close();
    }
    catch(std::ofstream::failure e){
        // fail to write then remove the key in pileup_data_files
        this->pileup_data_files.erase(chrom);
    }
    }
    
    return ;
}


pos_tc_t CallerFromAlignments::__pileup_treat_ctrl_a_chromosome_multiThread ( std::string chrom , bool uniqOnly )
{
    /*"""After this function is called, chr_pos_treat_ctrl will
     be reset and assigned to the pileup values of the given
     chromosome.
     
     """*/
    
    std::pair<std::vector<int>, std::vector<double> > treat_pv, ctrl_pv;
    pos_tc_t res_pos_tc;

    //long i;
    //float t;
    std::ifstream f;
    
    if (std::find(canoChromlist.begin(),canoChromlist.end(),chrom) == canoChromlist.end()) {
        error("chromosome is not valid : " + chrom);
        std::exit(1);
    }
    //assert chrom in self.chromosomes, "chromosome %s is not valid." % chrom;
    
    // check backup file of pileup values. If not exists, create
    // it. Otherwise, load them instead of calculating new pileup values.
/*    if (uniqOnly) {
        
        if (pileup_data_files.find( chrom ) != pileup_data_files.end())
        {
            
            try{
                f.open( pileup_data_files[ chrom ],std::ifstream::in );

                
                while(! f.eof()){
                    std::string line, ss_pos,ss_treat,ss_ctrl;
                    std::stringstream ss;
                    
                    if (! std::getline(f,line)){
                        break;
                    }
                    
                    ss << line;
                    std::getline(ss,ss_pos,'\t');
                    std::getline(ss,ss_treat,'\t');
                    std::getline(ss,ss_ctrl,'\t');
                    
                    res_pos_tc.pos.push_back(std::stol(ss_pos));
                    res_pos_tc.treat_v.push_back(std::stod(ss_treat));
                    res_pos_tc.ctrl_v.push_back(std::stod(ss_ctrl));
                    
                }
                
                f.close();
                
                
                return res_pos_tc;
            }
            catch(std::ifstream::failure e ){
                struct timeval start;
                gettimeofday(&start, NULL);
                std::string temp_filename = lexical_cast<std::string>(start.tv_sec) + "." + lexical_cast<std::string> ( start.tv_usec);
                //temp_fd, temp_filename = mkstemp();
                //os.close(temp_fd);
                pileup_data_files.insert(std::pair<std::string,std::string>(chrom,temp_filename));
            }
        }
        else{
            //temp_fd, temp_filename = mkstemp();
            //os.close(temp_fd);
            //pileup_data_files[ chrom ] = temp_filename;
            struct timeval start;
            gettimeofday(&start, NULL);
            std::string temp_filename = lexical_cast<std::string>(start.tv_sec) + "." + lexical_cast<std::string> ( start.tv_usec);
            pileup_data_files.insert(std::pair<std::string,std::string>(chrom,temp_filename));
            
        }
    }
    */
    
    std::vector<double> scale_factors;
    std::vector<int> ds;
    ds.push_back(this->d);
    scale_factors.push_back(treat_scaling_factor);
    
    if (PE_mode){
        
        treat_pv = treat->pileup_a_chromosome_pe ( chrom, scale_factors, 0.0 );
    }
    else{
        
        treat_pv = treat->pileup_a_chromosome( chrom, ds, scale_factors,uniqOnly, 0.0,true, this->end_shift );
        
        
    }
    if (not no_lambda_flag){
        if (PE_mode){
            // note, we pileup up PE control as SE control because
            // we assume the bias only can be captured at the
            // surrounding regions of cutting sites from control experiments.
            ctrl_pv = ctrl->pileup_a_chromosome_c( chrom, this->ctrl_d_s, this->ctrl_scaling_factor_s, this->lambda_bg );
        }
        else{
            
            
            ctrl_pv = ctrl->pileup_a_chromosome( chrom, this->ctrl_d_s, this->ctrl_scaling_factor_s,uniqOnly,
                                                this->lambda_bg,
                                                false );
        }
    }
    else{
        std::vector<int> pos;
        pos.push_back(treat_pv.first[treat_pv.first.size()-1]);
        std::vector<double> val;
        val.push_back(this->lambda_bg);
        
        ctrl_pv = std::pair<std::vector<int>,std::vector<double> > (pos,val);
        //   ctrl_pv = [treat_pv.first[][0][-1:], np.array([self.lambda_bg,], dtype="float32")]; // set a global lambda
    }
    
    __chrom_pair_treat_ctrl_multiThread( treat_pv, ctrl_pv,res_pos_tc);
    
    // clean treat_pv and ctrl_pv
    treat_pv.first.clear(); // = [];
    treat_pv.second.clear();
    ctrl_pv.first.clear(); //  = [];
    ctrl_pv.second.clear();
    

/*    if (uniqOnly) {
        
        try{
            std::ofstream f(this->pileup_data_files[ chrom ],std::ofstream::out);
            //    cPickle.dump( chr_pos_treat_ctrl, f , protocol=2 );
            
            for (size_t i=0 ; i < res_pos_tc.pos.size(); i ++ ) {
                
                //double tr_v = round(chr_pos_treat_ctrl.treat_v[i] * 100000)/100000;
                //double ct_v = round(chr_pos_treat_ctrl.ctrl_v[i] * 100000)/100000;
                
                double tr_v = res_pos_tc.treat_v[i];
                double ct_v = res_pos_tc.ctrl_v[i];
                
                f  << res_pos_tc.pos[i] << "\t" << tr_v << "\t" << std::setprecision(10) << ct_v << "\n";
            }
            f.close();
        }
        catch(std::ofstream::failure e){
            // fail to write then remove the key in pileup_data_files
            this->pileup_data_files.erase(chrom);
        }
    }*/
    
    return res_pos_tc;
}


std::vector<double> CallerFromAlignments::__cal_pscore ( std::vector<double> array1, std::vector<double> array2 ,bool uniqOnly)
{
    std::vector<double> s;
    
    if (array1.size() != array2.size()) {
        error("Two arraies shouble have the same length.");
        std::exit(1);
    }
    
    for (size_t i=0 ; i< array1.size();i++)
    {
      if(uniqOnly){
          double r = get_pscore( int(array1[i]), array2[i] );
          //std::cout << "cal_pscore " << array1[i] <<"\t"<< array2[i]<<"\t"<<r << std::endl;
          s.push_back(r);
          //s.push_back( get_pscore( int(array1[i]), array2[i] ));
      }
      else {
          double r = get_pscore_v2( int(array1[i]), array2[i] );

          s.push_back(r);
     }
    }
    return s;
}



void CallerFromAlignments::__cal_pvalue_qvalue_table ()
{
       /* """After this function is called, pqtable is built. All
        chromosomes will be iterated. So it will take some time.
        
        """*/

    Float64HashTable pvalue_stat(1) ;

    info( "#4 Start to calculate pvalue stat..." );

    long N = 0;
    int nhcal = 0;
    std::vector<double> keys;
    
    //for(size_t i =0; i< chromosomes.size(); i++ )
    debug("canochromlist size = " + std::to_string(canoChromlist.size()));
    for(size_t  i = 0; i < canoChromlist.size() ; i++)
    {
        
        std::string chrom = canoChromlist[ i ];
        
        int pre_p = 0;

        __pileup_treat_ctrl_a_chromosome( chrom , true); //unique reads only
        
        for (size_t j=0; j< chr_pos_treat_ctrl.pos.size() ; j++)
        {
            double this_v = get_pscore( int(chr_pos_treat_ctrl.treat_v[j]), chr_pos_treat_ctrl.ctrl_v[j] );
            
 
            int this_l = chr_pos_treat_ctrl.pos[j] - pre_p;
            
            N += this_l;
            
            double this_v_round = roundTop(this_v,10);
           
            if (pvalue_stat.has_key( this_v_round )){
                pvalue_stat.inc_item(this_v_round, this_l);
            }
            else{
                pvalue_stat.set_item(this_v_round, this_l);
                keys.push_back(this_v_round);
            }
            pre_p = chr_pos_treat_ctrl.pos[j];
        }
        //nhcal += chr_pos_treat_ctrl.pos.size();
        
    }
    
    
    long k = 1;//                           # rank
    double f = -1.0 * log10(N);
    //double pre_v = -2147483647;
    //int pre_l = 0;
    double pre_q = 2147483647; //             # save the previous q-value

    info( "#4 New pqtable ..." );
    pqtable = new Float64HashTable(1);
    
    //sort p-values from small to big, unique values
    std::sort(keys.begin(), keys.end(), std::greater<double>());
    
    
    for (auto key : keys)
    {
        double v = key; //p-value
        
        long l = (long)pvalue_stat.get_item(key);
        
        double q = v + (log10(k) + f);
        
        q = std::max(0.0,std::min(pre_q,q)); //           # make q-score monotonic
        
        pqtable->set_item(v,q);
        
        //std::cout << "set_item key  " << v << "\t" << q << std::endl;
        
        pre_q = q;
        k +=l;
        nhcal += 1;
    }
    
    info( " access pq hash for " + std::to_string(nhcal) + " times.");
    
    info("*** size of pqtable " + std::to_string(pqtable->size() ));

}

void CallerFromAlignments::__chromlist_call_peak_using_certain_criteria(std::vector<std::string> chromlist,PeakIO *peaks, call_peak_options_t copts )
{
    for ( auto chrom : chromlist){
        // treat/control bedGraph will be saved if requested by user.
        
        __chrom_call_peak_using_certain_criteria ( peaks, chrom, copts , this->save_bedGraph );
    }
}

void CallerFromAlignments::__add_broadpeak (PeakIO *peaks, std::string chrom, PeakContent lvl2peak, std::vector<PeakContent> lvl1peakset)
{
       /* """Internal function to create broad peak.

        *Note* lvl1peakset/strong_regions might be empty
        """*/
        //print lvl2peak["start"], lvl2peak["end"], lvl2peak["score"]
    int start      = lvl2peak.start;
    int end        = lvl2peak.end;

    if ( lvl1peakset.size() == 0)
    {
            /*try:
            # will complement by adding 1bps start and end to this region
            # may change in the future if gappedPeak format was improved.*/
        peaks->add(chrom, start, end, -1, lvl2peak.score, lvl2peak.pileup,
                       lvl2peak.pscore, lvl2peak.fc,
                   lvl2peak.qscore );
        return;

    }
    //Toadd blocks

    peaks->add(chrom, start, end, -1, lvl2peak.score,  lvl2peak.pileup,
               lvl2peak.pscore,  lvl2peak.fc,
               lvl2peak.qscore );
    
}

void CallerFromAlignments::__chromlist_call_broadpeak_using_certain_criteria(std::vector<std::string> chromlist,PeakIO *peaks, call_peak_options_t copts )
{
    for ( auto chrom : chromlist){
        // treat/control bedGraph will be saved if requested by user.
        
        PeakIO *lvl1peaks = new PeakIO();
        PeakIO *lvl2peaks = new PeakIO();
        
        __chrom_call_broadpeak_using_certain_criteria ( lvl1peaks,lvl2peaks, chrom, copts , this->save_bedGraph  );
        
        
        
        
        //now combine lvl1 and lvl2 peaks
        std::vector<std::string> chrs = lvl1peaks->get_chr_names();
        //broadpeaks = BroadPeakIO();
        
        info("linking regions between lvl1_peaks");
        // use lvl2_peaks as linking regions between lvl1_peaks
        for(auto chrom : chrs)
        {
            //int tmp_n = 0;
            std::vector<PeakContent> lvl1peakschrom = lvl1peaks->get_data_from_chrom(chrom);
            std::vector<PeakContent> lvl2peakschrom = lvl2peaks->get_data_from_chrom(chrom);
            unsigned int j = 0;
            unsigned int i = 0;
            std::vector<PeakContent> tmppeakset;            // to temporarily store lvl1 region inside a lvl2 region
            // our assumption is lvl1 regions should be included in lvl2 regions
                for(;i < lvl2peakschrom.size() ; i ++ )
                {
                    // for each lvl2 peak, find all lvl1 peaks inside
                    // I assume lvl1 peaks can be ALL covered by lvl2 peaks.
                    PeakContent lvl2 = lvl2peakschrom[i];

                    while(j < lvl1peakschrom.size())
                    {
                        PeakContent lvl1 = lvl1peakschrom[j];
                        j += 1;
                        if ( lvl2.start <= lvl1.start && lvl1.end <= lvl2.end)
                        {
                            tmppeakset.push_back(lvl1);
                        }
                        else{
                            // make a hierarchical broad peak
                            //print lvl2["start"], lvl2["end"], lvl2["score"]
                            //info("add broad peak 1");
                            __add_broadpeak ( peaks, chrom, lvl2, tmppeakset);
                            
                            tmppeakset.clear();
                            break;
                        }
                    }//end while
                    if (j >= lvl1peakschrom.size())
                    {
                        // no more strong (aka lvl1) peaks left
                        //info("add broad peak 2");
                        __add_broadpeak ( peaks, chrom, lvl2, tmppeakset);
                        
                        tmppeakset.clear();
                        break;
                    }
                }//end for

            // add the rest lvl2 peaks
            for(unsigned int k = i + 1; k > lvl2peakschrom.size() ; k ++ )
            {
                //info("add broad peak 3");
                __add_broadpeak( peaks, chrom, lvl2peakschrom[k], tmppeakset );
            }
            
        }//end for chrom
        
        delete lvl1peaks;
        delete lvl2peaks;
    }// end all chroms
}

void CallerFromAlignments::merge_peaklist(PeakIO *peaks, std::vector<PeakIO *> peaks_list,bool uniqOnly)
{
    //if(uniqOnly)
    //{
        for (auto p : peaks_list ) {
            
            peaks->merge_peaks(p);
            
        }
   
}

void CallerFromAlignments::call_peaks ( PeakIO *peaks, call_peak_options_t copts )
{
        /*"""Call peaks for all chromosomes. Return a PeakIO object.
        
        """*/
    
    
    //std::string chrom;
    std::string s;
    //bytes tmp_bytes;

    //prepare p-q table
    if (pqtable == nullptr && copts.uniqOnly){
        info("#4 Pre-compute pvalue-qvalue table...");
        if (copts.auto_cutoff){
            info("#4 Cutoff will be automatically decided!");
            //ADD __pre_computes( max_gap, min_length );
        }
        else{
            __cal_pvalue_qvalue_table();
        }
    }
    
    info("#4 Call peaks for each chromosome...");
    
    //multi-threading
    std::vector<PeakIO *> peaks_list;
    
    //each thread is assigned a PeakIO object
    for (int i=0; i < copts.threads_to_use; i++)
    {
        peaks_list.push_back(new PeakIO());
    }
    
    std::vector<std::thread *> t;
    
    //int total_num_chroms = chromosomes.size();
    int total_num_chroms = canoChromlist.size();
    
    info("total number of chroms " + std::to_string(total_num_chroms));
    int chrom_num_per_thread = total_num_chroms/copts.threads_to_use;

    //auto start = std::chrono::high_resolution_clock::now();
    int leftOver_num_chroms = total_num_chroms - copts.threads_to_use * chrom_num_per_thread;
    
    for (int k =0; k < leftOver_num_chroms; k ++) //assign chroms to the first threads with (chrom_num_per_thread + 1) chroms
    {
        std::vector<std::string> chromlist;
        int num_chroms =chrom_num_per_thread + 1;
        
        info("num_chroms per thread " + std::to_string(num_chroms));
        
        for (int i= 0; i < num_chroms; i++) {
            chromlist.push_back(canoChromlist[i + k * num_chroms]);
        }
        
        t.push_back(new std::thread(&CallerFromAlignments::__chromlist_call_peak_using_certain_criteria,this,chromlist,peaks_list[k],copts));
    }
    

    for (int k = leftOver_num_chroms; k < copts.threads_to_use; k++) {
        std::vector<std::string> chromlist;
        int num_chroms =chrom_num_per_thread;
        
        info("num_chroms per thread " + std::to_string(num_chroms));
        
        for (int i= 0; i < num_chroms; i++) {
            chromlist.push_back(canoChromlist[i + leftOver_num_chroms + k*chrom_num_per_thread]);
        }
        
        t.push_back(new std::thread(&CallerFromAlignments::__chromlist_call_peak_using_certain_criteria,this,chromlist,peaks_list[k],copts));
        
        if (num_chroms != chrom_num_per_thread) {
            break;
        }
    }
    info("before threads join");
    for (int k =0; k < copts.threads_to_use; k++) {
        t[k]->join();
    }
    info("after threads join");
    //merge results
    merge_peaklist(peaks,peaks_list,copts.uniqOnly);
    
    info("after merge peaks");
    //clean up
    
    for (int i = 0; i < copts.threads_to_use; i++)
    {
        delete t[i];
        delete peaks_list[i];
    }
    
}

void CallerFromAlignments::call_broadpeaks ( PeakIO *peaks, call_peak_options_t copts )
{
    
    /*"""Call peaks for all chromosomes. Return a PeakIO object."""*/
    
    //std::string chrom;
    std::string s;
    //bytes tmp_bytes;
    
    //prepare p-q table
    if (pqtable == nullptr && copts.uniqOnly){
        info("#4 Pre-compute pvalue-qvalue table...");
        if (copts.auto_cutoff){
            info("#4 Cutoff will be automatically decided!");
            //ADD __pre_computes( max_gap, min_length );
        }
        else{
            __cal_pvalue_qvalue_table();
        }
    }
    
    info("#4 Call peaks for each chromosome...");
    
    //multi-threading
    std::vector<PeakIO *> peaks_list;
    //std::vector<PeakIO *> peaks_list2;
    
    //each thread is assigned a PeakIO object
    for (int i=0; i < copts.threads_to_use; i++)
    {
        peaks_list.push_back(new PeakIO());
    }
    
    std::vector<std::thread *> t;
    
    //int total_num_chroms = chromosomes.size();
    int total_num_chroms = canoChromlist.size();
    
    info("total number of chroms " + std::to_string(total_num_chroms));
    int chrom_num_per_thread = total_num_chroms/copts.threads_to_use;
    
    //auto start = std::chrono::high_resolution_clock::now();
    int leftOver_num_chroms = total_num_chroms - copts.threads_to_use * chrom_num_per_thread;
    
    for (int k =0; k < leftOver_num_chroms; k ++) //assign chroms to the first threads with (chrom_num_per_thread + 1) chroms
    {
        std::vector<std::string> chromlist;
        int num_chroms =chrom_num_per_thread + 1;
        
        info("num_chroms per thread " + std::to_string(num_chroms));
        
        for (int i= 0; i < num_chroms; i++) {
            chromlist.push_back(canoChromlist[i + k * num_chroms]);
        }
        
        t.push_back(new std::thread(&CallerFromAlignments::__chromlist_call_broadpeak_using_certain_criteria,this,chromlist,peaks_list[k],copts));
    }
    
    
    for (int k = leftOver_num_chroms; k < copts.threads_to_use; k++) {
        std::vector<std::string> chromlist;
        int num_chroms =chrom_num_per_thread;
        
        info("num_chroms per thread " + std::to_string(num_chroms));
        
        for (int i= 0; i < num_chroms; i++) {
            chromlist.push_back(canoChromlist[i + leftOver_num_chroms + k*chrom_num_per_thread]);
        }
        
        t.push_back(new std::thread(&CallerFromAlignments::__chromlist_call_broadpeak_using_certain_criteria,this,chromlist,peaks_list[k],copts));
        
        if (num_chroms != chrom_num_per_thread) {
            break;
        }
    }
    info("before threads join");
    for (int k =0; k < copts.threads_to_use; k++) {
        t[k]->join();
    }
    info("after threads join");
    //merge results
    merge_peaklist(peaks,peaks_list,copts.uniqOnly);
    
    info("after merge peaks");
    //clean up
    

    for (int i = 0; i < copts.threads_to_use; i++)
    {
        delete t[i];
        delete peaks_list[i];
    }
    
}


bool CallerFromAlignments::__close_peak_wo_subpeaks (bool uniqOnly,std::vector<peak_content_t> peak_content, PeakIO *peaks, int min_length, std::string chrom, int smoothlen, std::vector<std::vector<double> > score_array_s, std::vector<double> score_cutoff_s)
{
       /* """Close the peak region, output peak boundaries, peak summit
        and scores, then add the peak to peakIO object. peak_content contains [start, end, treat_p, ctrl_p, index_in_score_array], peaks: a PeakIO object
        """*/

    int peak_length = peak_content[peak_content.size()-1].e - peak_content[0].s;
    
    if (peak_length >= min_length)
    { // # if the peak is too small, reject it
        std::vector<int> tsummit;
        std::vector<int> tsummit_index;
        int summit_pos   = 0;
        double summit_value = 0.0;
        double summit_p_score = 0.0;
        double summit_q_score = 0.0;
        double summit_treat = peak_content[ int((peak_content.size() + 1) / 2) - 1 ].t_val;
        double summit_ctrl = peak_content[ int((peak_content.size() + 1) / 2) - 1 ].c_val;
        
        if (not uniqOnly) {
            //for candidate multi-mapper peak regions, we just output the range
            peaks->add( chrom,           // chromosome
                       peak_content[0].s, // start
                       peak_content[peak_content.size()-1].e, // end
                       summit_pos, // summit position
                       summit_q_score, // score at summit
                       summit_treat, // pileup
                       summit_p_score, // pvalue
                       1.0* ( summit_treat + pseudocount ) / (summit_ctrl + pseudocount ), // fold change
                       summit_q_score // qvalue
                       );
            return true;
            
        }
         
        for(size_t i=0; i< peak_content.size();i++)
        {
            int tstart = peak_content[i].s;
            int tend = peak_content[i].e;
            double ttreat_p = peak_content[i].t_val;
           // double tctrl_p = peak_content[i].c_val;
            //int tlist_scores_p = peak_content[i].index;
            
            double tscore = ttreat_p; // # use qscore as general score to find summit
            if(summit_value == 0 or summit_value < tscore){
                tsummit.clear();
                tsummit_index.clear();
                tsummit.push_back((tend + tstart) / 2);
                tsummit_index.push_back(i);
                summit_value = tscore;
            }
            else {
                if (summit_value == tscore){
                    //# remember continuous summit values
                    tsummit.push_back(int((tend + tstart) / 2));
                    tsummit_index.push_back( i );
                }
            }
        }
            // the middle of all highest points in peak region is defined as summit
        int midindex = int((tsummit.size() + 1) / 2) - 1;
        summit_pos    = tsummit[ midindex ];
        int summit_index  = tsummit_index[ midindex ];
        
        summit_treat = peak_content[ summit_index ].t_val;
        summit_ctrl = peak_content[ summit_index ].c_val;
        
        //# this is a double-check to see if the summit can pass cutoff values.
        for (size_t i=0 ; i < score_cutoff_s.size();i++)
        {
            if (score_cutoff_s[i] > score_array_s[ i ][ peak_content[ summit_index ].index ])
            {
                return false; // # not passed, then disgard this peak.
            }
        }
        
           summit_p_score = get_pscore( int(summit_treat), summit_ctrl );
            double summit_p_score_round = roundTop(summit_p_score,10);
            summit_q_score = pqtable->get_item(summit_p_score_round);

            peaks->add( chrom,           // chromosome
                       peak_content[0].s, // start
                       peak_content[peak_content.size()-1].e, // end
                       summit_pos, // summit position
                       summit_q_score, // score at summit
                       summit_treat, // pileup
                       summit_p_score, // pvalue
                       1.0* ( summit_treat + pseudocount ) / (summit_ctrl + pseudocount ), // fold change
                       summit_q_score // qvalue
                      );

        
        //# start a new peak
        return true;
    }
    else {
        return false;
    }
}

bool CallerFromAlignments::__close_peak_with_subpeaks (std::vector<peak_content_t> peak_content, PeakIO *peaks, int min_length, std::string chrom, int smoothlen, std::vector<std::vector<double> > score_array_s, std::vector<double> score_cutoff_s,float min_valley )
{
       /* """Algorithm implemented by Ben, to profile the pileup signals
        within a peak region then find subpeak summits. This method is
        highly recommended for TFBS or DNAase I sites.
        
        """*/

    info("__close_peak_with_subpeaks");
    int peak_length = peak_content[ peak_content.size()-1].e - peak_content[ 0 ].s;
    
    if (peak_length < min_length){
        return false;
    }// if the region is too small, reject it

    //Add 10 bp padding to peak region so that we can get true minima
    int end = peak_content[peak_content.size()-1].e + 10;
    int start = peak_content[0].s - 10;
    int start_boundary = 10; //// # this is the offset of original peak boundary in peakdata list.

    
    if (start < 0){
        start_boundary = 10 + start;// # this is the offset of original peak boundary in peakdata list.
        start = 0;
    }

    std::vector<double> peakdata(end-start,0); // # save the scores (qscore) for each position in this region
    std::vector<int> peakindices(end-start,0); // # save the indices for each position in this region
    for (size_t i=0;i < peak_content.size();i++)
    {
        //(tstart, tend, ttreat_p, tctrl_p, tlist_scores_p) = peak_content[i];
        //tscore = self.pqtable[ get_pscore(int(ttreat_p), tctrl_p) ] # use qscore as general score to find summit
        int tscore = peak_content[i].t_val;// # use pileup as general score to find summit
        int m = peak_content[i].s - start + start_boundary;
        int n = peak_content[i].e - start + start_boundary;
        for (int k = m; k< n; k++) {
            peakdata[k] = tscore;
            peakindices[k] = i;
        }
    }
    std::vector<int> summit_offsets;
    maxima(peakdata,summit_offsets,smoothlen); // # offsets are the indices for summits in peakdata/peakindices array.
        //print "maxima:",summit_offsets
    if (summit_offsets.size() == 0){
            // **failsafe** if no summits, fall back on old approach #
        return __close_peak_wo_subpeaks(true,peak_content, peaks, min_length, chrom, smoothlen, score_array_s, score_cutoff_s);
    }
    else{
        // remove maxima that occurred in padding
        //m = np.searchsorted(summit_offsets, start_boundary);
        std::sort(summit_offsets.begin(),summit_offsets.end());
        std::vector<int>::iterator m = std::lower_bound(summit_offsets.begin(),summit_offsets.end(),start_boundary);
        std::vector<int>::iterator n = std::upper_bound(summit_offsets.begin(),summit_offsets.end(), peak_length + start_boundary);
        
        
        summit_offsets.erase (n+1,summit_offsets.end());
        summit_offsets.erase(summit_offsets.begin(),m-1);
        //summit_offsets = summit_offsets[m:n];
    }
    
    summit_offsets = enforce_peakyness(peakdata, summit_offsets);
    
        //print "enforced:",summit_offsets
    if (summit_offsets.size() == 0){
            // **failsafe** if no summits, fall back on old approach #
        return __close_peak_wo_subpeaks(true,peak_content, peaks, min_length, chrom, smoothlen, score_array_s, score_cutoff_s);
    }
    
    for (size_t k=0; k < summit_offsets.size(); k++){

        int summit_index = peakindices[summit_offsets[k]];
        int summit_offset = summit_offsets[k] - start_boundary;

        double summit_treat = peak_content[ summit_index ].t_val;
        double summit_ctrl = peak_content[ summit_index ].c_val;

        double summit_p_score = get_pscore( int(summit_treat), summit_ctrl );
        double summit_p_score_round = roundTop(summit_p_score,10);
        double summit_q_score = pqtable->get_item(summit_p_score_round);
        
        
        for(size_t i=0; i < score_cutoff_s.size();i++)
        {
            if (score_cutoff_s[i] > score_array_s[ i ][ peak_content[ summit_index ].index ]){
                return false; // # not passed, then disgard this summit.
            }
        }

        peaks->add( chrom,
                       peak_content[ 0 ].s,
                       peak_content[peak_content.size()-1].e,
                       start + summit_offset,
                       summit_q_score,
                       summit_treat,
                       summit_p_score,
                       1.0* ( summit_treat + pseudocount ) / ( summit_ctrl + pseudocount ), //# fold change
                       summit_q_score
                      );
    }
    //start a new peak
    return true;
}

std::vector<double> CallerFromAlignments::__cal_qscore ( std::vector<double> array1, std::vector<double> array2 )
{
    std::vector<double> s;

    if(array1.size() != array2.size())
    {
        error("tow lists should have the same length.");
        std::exit(1);
    }
    
    for (size_t i=0; i< array1.size(); i++)
    {
        //std::cout.precision(10);
        //std::cout << int(array1[i]) << std::fixed << "\t" << array2[i] << std::endl;
        
        double p = get_pscore( int(array1[i]), array2[i] );
        
        //std::cout << "p = " << std::fixed << p << std::endl;
        double p_round = roundTop(p,10);
        
        //std::cout << "p = " << std::fixed << p_round << std::endl;
        
        double q = pqtable->get_item(p_round);
        
        
        s.push_back(q);

    }
   
    return s;
}

void CallerFromAlignments::__write_bedGraph_for_a_chromosome(std::string chrom,pos_tc_t chr_pos_treat_ctrl_multiThread)
{
    //"""Write treat/control values for a certain chromosome into a
    //    specified file handler.
    //    """
    std::vector<int> pos_array;
    std::vector<double> treat_array;
    std::vector<double> ctrl_array;
    int l;
    int p,pre_p_t,pre_p_c; //current position, previous position for treat, previous position for control
    double pre_v_t, pre_v_c, v_t, v_c; // previous value for treat, for control, current value for treat and control
    
    if(chrom.length() > 6)
    {
        return;
    }
    //double denominator = 1; //
    
    pos_array = chr_pos_treat_ctrl_multiThread.pos;
    treat_array = chr_pos_treat_ctrl_multiThread.treat_v;
    ctrl_array = chr_pos_treat_ctrl_multiThread.ctrl_v;
    
    std::ofstream tmp_bedGraph_treat_f;
    std::ofstream tmp_bedGraph_ctrl_f;
    
    tmp_bedGraph_treat_f.open(this->bedGraph_filename_prefix + "_treat_uniq_"+chrom+".bdg",std::ofstream::out);
    
    tmp_bedGraph_ctrl_f.open(this->bedGraph_filename_prefix + "_ctrl_uniq_"+chrom+".bdg",std::ofstream::out);
    
    tmp_bedGraph_treat_f << "track type=bedGraph name=\"treatment pileup\" description=\"treatment pileup after possible scaling for \'UTF-8\'" << "\n";
    
    tmp_bedGraph_ctrl_f << "track type=bedGraph name=\"control lambda\" description=\"control lambda after possible scaling for \'UTF-8\'" << "\n";
    
    
    l = pos_array.size();
    if(l ==0)
    {
        return;
    }
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
            if(pre_v_t > 0)
            {
            tmp_bedGraph_treat_f << chrom << "\t" << pre_p_t << "\t" << p << "\t" << std::setprecision(5) << pre_v_t << "\n";
            }
            pre_v_t = v_t;
            pre_p_t = p;
        }
        if ( std::abs(pre_v_c - v_c) > 1e-5) // precision is 5 digits
        {
            if(pre_v_c > 0)
            {
            tmp_bedGraph_ctrl_f << chrom << "\t" << pre_p_c << "\t" << p << "\t" << std::setprecision(5) << pre_v_c << "\n";
            }
            pre_v_c = v_c;
            pre_p_c = p;
        }
            
    }
    p = pos_array[l-1];
    if(pre_v_t > 0) {
        tmp_bedGraph_treat_f << chrom << "\t" << pre_p_t << "\t" << p << "\t" << std::setprecision(5) << pre_v_t << "\n";
    }
    if(pre_v_c > 0)
    {
        tmp_bedGraph_ctrl_f << chrom << "\t" << pre_p_c << "\t" << p << "\t" << std::setprecision(5) << pre_v_c << "\n";
    }
    
    tmp_bedGraph_treat_f.close();
    tmp_bedGraph_ctrl_f.close();

}
void CallerFromAlignments::__chrom_call_peak_using_certain_criteria ( PeakIO *peaks, std::string chrom, call_peak_options_t copts,  bool save_bedGraph )
{
        /*""" Call peaks for a chromosome.

        Combination of criteria is allowed here.

        peaks: a PeakIO object, the return value of this function
        scoring_function_s: symbols of functions to calculate score as score=f(x, y) where x is treatment pileup, and y is control pileup
        save_bedGraph     : whether or not to save pileup and control into a bedGraph file
        """*/

    //list score_array_s ;         //# list to keep different types of scores
    std::vector<peak_content_t> peak_content_list  ;         //#  to store information for a
                                       /* #  chunk in a peak region, it
                                        #  contains lists of: 1. left
                                        #  position; 2. right
                                        #  position; 3. treatment
                                        #  value; 4. control value;
                                        #  5. list of scores at this
                                        #  chunk*/
    //long  lastp;
    //float tp, cp;
    
    if (copts.scoring_function_symbols.size() != copts.score_cutoff_s.size()) {
        error("number of functions and cutoffs should be the same!");
        std::exit(1);
    }
    //peak_content = [] ; //          # to store points above cutoff

    //# first, build pileup, chr_pos_treat_ctrl
    //# this step will be speeped up if pqtable is pre-computed.
    
    pos_tc_t chr_pos_treat_ctrl_multiThread = __pileup_treat_ctrl_a_chromosome_multiThread( chrom,copts.uniqOnly );

    //# while save_bedGraph is true, invoke __write_bedGraph_for_a_chromosome
    if (save_bedGraph){
        __write_bedGraph_for_a_chromosome ( chrom,chr_pos_treat_ctrl_multiThread );
        //info("write bedgraph for a chromosome.");
    }
    
    //keep all types of scores needed
    info("#4 after pileup ... " + chrom);
    //t0 = ttime()
    std::vector<std::vector<double> > score_array_s ;
    for( size_t i=0;i< copts.scoring_function_symbols.size();i++)
    {
        if (copts.scoring_function_symbols[i] == "p"){
            score_array_s.push_back(__cal_pscore( chr_pos_treat_ctrl_multiThread.treat_v, chr_pos_treat_ctrl_multiThread.ctrl_v, copts.uniqOnly ));
        }
        if(copts.scoring_function_symbols[i] == "q")
        {
            score_array_s.push_back(__cal_qscore( chr_pos_treat_ctrl_multiThread.treat_v, chr_pos_treat_ctrl_multiThread.ctrl_v )) ;
        }
    }

    std::vector<int> above_cutoff_startpos, above_cutoff_endpos;
    std::vector<int> above_cutoff_index_array;
    std::vector<int> FE_ret;

    apply_multiple_cutoffs(score_array_s,copts.score_cutoff_s,above_cutoff_index_array);

    
    //info("#4 after apply multiple cutoffs ..." + std::to_string( above_cutoff_index_array.size() ) );
    if (above_cutoff_index_array.size() == 0){
            //# nothing above cutoff
        return ;
    }
    
    for (size_t j =0; j < above_cutoff_index_array.size(); j++) {
    above_cutoff_endpos.push_back(chr_pos_treat_ctrl_multiThread.pos[above_cutoff_index_array[j]]);
    above_cutoff_startpos.push_back(chr_pos_treat_ctrl_multiThread.pos[above_cutoff_index_array[j]-1]);
        
    }
    if (above_cutoff_index_array[0] == 0){
        //# first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
        above_cutoff_startpos[0] = 0;
    }

    peak_content_t pt ;
    
    pt.s =above_cutoff_startpos[0];
    pt.e =above_cutoff_endpos[0];
    pt.index = above_cutoff_index_array[0];
    pt.t_val = chr_pos_treat_ctrl_multiThread.treat_v[pt.index];
    pt.c_val = chr_pos_treat_ctrl_multiThread.ctrl_v[pt.index];
    
    peak_content_list.push_back(pt);
    
    long lastp = pt.e;
    //info("#4 above_cutoff_startpos.size " + std::to_string(above_cutoff_startpos.size()));
    for (size_t i =1; i < above_cutoff_startpos.size(); i++ )
    {
        peak_content_t pt2;
        
        pt2.s = above_cutoff_startpos[ i ];
        pt2.e = above_cutoff_endpos[ i ];
        pt2.index = above_cutoff_index_array[ i ];
       // acs_ptr += 1;
        //ace_ptr += 1;
        //acia_ptr+= 1;
        pt2.t_val = chr_pos_treat_ctrl_multiThread.treat_v[ pt2.index ];
        pt2.c_val = chr_pos_treat_ctrl_multiThread.ctrl_v[ pt2.index ];
        long tl = pt2.s - lastp;
        
        if (tl <= copts.max_gap){
                // append.
            peak_content_list.push_back(pt2 );
            lastp = pt2.e; // #above_cutoff_endpos[i]
        }
        else{
           //info("#4 before close peaks 1...");
                //close
            if (copts.call_summits){
                __close_peak_with_subpeaks (peak_content_list, peaks, copts.min_length, chrom, copts.min_length, score_array_s, copts.score_cutoff_s ); // # smooth length is min_length, i.e. fragment size 'd'
            }
            else{
                //info("#4 before close peaks 2 ...");
                __close_peak_wo_subpeaks (copts.uniqOnly,peak_content_list, peaks, copts.min_length, chrom, copts.min_length, score_array_s, copts.score_cutoff_s ); //# smooth length is min_length, i.e. fragment size 'd'
                //info("#4 after close peaks 2 ... ");
             }//end else
            peak_content_list.clear();
            peak_content_list.push_back(pt2);
            lastp = pt2.e;// #above_cutoff_endpos[i];
        }//end else
    } //end for
    
    //save the last peak
    //info(" peaks size " + std::to_string(peaks->peaks[chrom].size()));
         
    if (peak_content_list.size()==0){
        return ;
    }
    else{
        if (copts.call_summits){
            __close_peak_with_subpeaks (peak_content_list, peaks, copts.min_length, chrom, copts.min_length, score_array_s, copts.score_cutoff_s ); // # smooth length is min_length, i.e. fragment size 'd'
        }
        else{
            __close_peak_wo_subpeaks (copts.uniqOnly,peak_content_list, peaks, copts.min_length, chrom, copts.min_length, score_array_s, copts.score_cutoff_s ); //# smooth length is min_length, i.e. fragment size 'd'

        }
    }
    //info("close peaks chrom = " + chrom);
    //#print "close peaks -- chrom:",chrom,"  time:", ttime() - t0
   // return peaks;
}

//call broad peak
void CallerFromAlignments::__chrom_call_broadpeak_using_certain_criteria ( PeakIO *peaks, PeakIO *peaks2,std::string chrom, call_peak_options_t copts,  bool save_bedGraph )
{
    /*""" Call peaks for a chromosome.
     
     Combination of criteria is allowed here.
     
     peaks: a PeakIO object, the return value of this function
     scoring_function_s: symbols of functions to calculate score as score=f(x, y) where x is treatment pileup, and y is control pileup
     save_bedGraph     : whether or not to save pileup and control into a bedGraph file
     """*/
    
    //list score_array_s ;         //# list to keep different types of scores
    std::vector<peak_content_t> peak_content_list  ;         //#  to store information for a
    
    if (copts.scoring_function_symbols.size() != copts.score_cutoff_s.size()) {
        error("number of functions and cutoffs should be the same!");
        std::exit(1);
    }
    if(copts.lvl1_cutoff_s.size() != copts.lvl2_cutoff_s.size())
    {
        error("number of functions and cutoffs should be the same!");
        std::exit(1);
    }
    
    //# first, build pileup, chr_pos_treat_ctrl
    //# this step will be speeped up if pqtable is pre-computed.
    
    pos_tc_t chr_pos_treat_ctrl_multiThread = __pileup_treat_ctrl_a_chromosome_multiThread( chrom,copts.uniqOnly );
    
    //# while save_bedGraph is true, invoke __write_bedGraph_for_a_chromosome
    if (save_bedGraph){
        __write_bedGraph_for_a_chromosome ( chrom,chr_pos_treat_ctrl_multiThread );
        //info("write bedgraph for a chromosome.");
    }
    
    //keep all types of scores needed
    std::vector<std::vector<double> > score_array_s ;
    for( size_t i=0;i< copts.scoring_function_symbols.size();i++)
    {
        if (copts.scoring_function_symbols[i] == "p"){
            score_array_s.push_back(__cal_pscore( chr_pos_treat_ctrl_multiThread.treat_v, chr_pos_treat_ctrl_multiThread.ctrl_v, copts.uniqOnly ));
        }
        if(copts.scoring_function_symbols[i] == "q")
        {
            score_array_s.push_back(__cal_qscore( chr_pos_treat_ctrl_multiThread.treat_v, chr_pos_treat_ctrl_multiThread.ctrl_v )) ;
        }
    }
    
    std::vector<int> above_cutoff_startpos, above_cutoff_endpos;
    std::vector<int> above_cutoff_index_array;
    std::vector<int> FE_ret;
 

    apply_multiple_cutoffs(score_array_s,copts.lvl1_cutoff_s,above_cutoff_index_array);

    
    info("#4 after apply multiple cutoffs ..."  );
    if (above_cutoff_index_array.size() == 0){
        //# nothing above cutoff
        return ;
    }
    
    for (size_t j =0; j < above_cutoff_index_array.size(); j++) {
        above_cutoff_endpos.push_back(chr_pos_treat_ctrl_multiThread.pos[above_cutoff_index_array[j]]);
        above_cutoff_startpos.push_back(chr_pos_treat_ctrl_multiThread.pos[above_cutoff_index_array[j]-1]);
        
    }
    if (above_cutoff_index_array[0] == 0){
        //# first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
        above_cutoff_startpos[0] = 0;
    }
    
    peak_content_t pt ;
    
    pt.s =above_cutoff_startpos[0];
    pt.e =above_cutoff_endpos[0];
    pt.index = above_cutoff_index_array[0];
    pt.t_val = chr_pos_treat_ctrl_multiThread.treat_v[pt.index];
    pt.c_val = chr_pos_treat_ctrl_multiThread.ctrl_v[pt.index];
    
    peak_content_list.push_back(pt);
    
    int lastp = pt.e;
    //info("#4 above_cutoff_startpos.size " + std::to_string(above_cutoff_startpos.size()));
    for (size_t i =1; i < above_cutoff_startpos.size(); i++ )
    {
        peak_content_t pt2;
        
        pt2.s = above_cutoff_startpos[ i ];
        pt2.e = above_cutoff_endpos[ i ];
        pt2.index = above_cutoff_index_array[ i ];
        // acs_ptr += 1;
        //ace_ptr += 1;
        //acia_ptr+= 1;
        pt2.t_val = chr_pos_treat_ctrl_multiThread.treat_v[ pt2.index ];
        pt2.c_val = chr_pos_treat_ctrl_multiThread.ctrl_v[ pt2.index ];
        long tl = pt2.s - lastp;
        
        if (tl <= copts.lvl1_max_gap){
            // append.
            peak_content_list.push_back(pt2 );
            lastp = pt2.e; // #above_cutoff_endpos[i]
            
        }
        else{
            //info("#4 before close peaks 1...");
            __close_peak_for_broad_region (copts.uniqOnly,peak_content_list, peaks, copts.min_length, chrom, copts.min_length, score_array_s ); // # smooth
            peak_content_list.clear();
            peak_content_list.push_back(pt2);
            lastp = pt2.e;// #above_cutoff_endpos[i];
            
        }//end else
    } //end for
    
    if (peak_content_list.size()==0){
        return ;
    }
    else{
        //info("close peaks 2 ...");
        __close_peak_for_broad_region (copts.uniqOnly,peak_content_list, peaks, copts.min_length, chrom, copts.min_length, score_array_s ); // # smooth length is min_length,
    }
    
    //*********************peaks2: weak peaks
    //lvl2 : weak peaks
    info("=======start weak peaks =====================");
    
    peak_content_list.clear(); // to store points above cutoff
    above_cutoff_index_array.clear();
    above_cutoff_endpos.clear();
    above_cutoff_startpos.clear();
    
    //get the regions with scores above cutoffs
    apply_multiple_cutoffs(score_array_s,copts.lvl2_cutoff_s,above_cutoff_index_array);
    
    if (above_cutoff_index_array.size() == 0){
        //# nothing above cutoff
        return ;
    }
    for (size_t j =0; j < above_cutoff_index_array.size(); j++) {
        above_cutoff_endpos.push_back(chr_pos_treat_ctrl_multiThread.pos[above_cutoff_index_array[j]]);
        above_cutoff_startpos.push_back(chr_pos_treat_ctrl_multiThread.pos[above_cutoff_index_array[j]-1]);
        
    }
    if (above_cutoff_index_array[0] == 0){
        //# first element > cutoff, fix the first point as 0. otherwise it would be the last item in data[chrom]['pos']
        above_cutoff_startpos[0] = 0;
    }
 
    //first bit of region above cutoff
    //peak_content_t pt ;
    
    pt.s =above_cutoff_startpos[0];
    pt.e =above_cutoff_endpos[0];
    pt.index = above_cutoff_index_array[0];
    pt.t_val = chr_pos_treat_ctrl_multiThread.treat_v[pt.index];
    pt.c_val = chr_pos_treat_ctrl_multiThread.ctrl_v[pt.index];
    
    peak_content_list.push_back(pt);

    lastp = pt.e;

    for (size_t i =1; i < above_cutoff_startpos.size(); i++ )
    {
        // for everything above cutoff
        peak_content_t pt2;
        
        pt2.s = above_cutoff_startpos[ i ];
        pt2.e = above_cutoff_endpos[ i ];
        pt2.index = above_cutoff_index_array[ i ];
        // acs_ptr += 1;
        //ace_ptr += 1;
        //acia_ptr+= 1;
        pt2.t_val = chr_pos_treat_ctrl_multiThread.treat_v[ pt2.index ];
        pt2.c_val = chr_pos_treat_ctrl_multiThread.ctrl_v[ pt2.index ];
        long tl = pt2.s - lastp;
        
            if (tl <= copts.lvl2_max_gap)
            {
                // append
                peak_content_list.push_back(pt2 );
                lastp = pt2.e;
            }
            else{
                // close
                 __close_peak_for_broad_region (copts.uniqOnly,peak_content_list, peaks2, copts.min_length, chrom, copts.lvl2_max_gap/2, score_array_s );
                //info("close peak 1");
                
                peak_content_list.clear();
                peak_content_list.push_back(pt2);
                lastp = pt2.e;
            }
        }
        // save the last peak
        if (peak_content_list.size()==0){
                return ;
        }
        else{
                __close_peak_for_broad_region (copts.uniqOnly,peak_content_list, peaks2, copts.min_length, chrom, copts.lvl2_max_gap/2, score_array_s ); // # smooth length is min_length,
            //info("close peak 2");
        }
}

std::vector<double> CallerFromAlignments::__cal_FE ( std::vector<double> array1, std::vector<double> array2 )
{
    std::vector<double> s;
    
    if (array1.size() != array2.size()) {
        error("Two arraies shouble have the same length.");
        std::exit(1);
    }
    
    for (size_t i=0 ; i< array1.size();i++)
    {
        double fc = (array1[0]+1.0)/(array2[0]+1.0);
        s.push_back(fc);
        
    }
    return s;
}

double CallerFromAlignments::mean_from_value_length(std::vector<double> array1, std::vector<int> lenlist)
{
    int total_len =0;
    double sum = 0;
    
    if(array1.size() != lenlist.size())
    {
        error("value list and length list shouble have the same length.");
        std::exit(1);
    }
    if (array1.size() == 0) {
        return 1.0;
    }
    for (unsigned int i=0; i < lenlist.size(); i++) {
        total_len += lenlist[i];
        sum += 1.0 * array1[i] * lenlist[i];
    }
    
    return sum/total_len;
}
bool CallerFromAlignments::__close_peak_for_broad_region (bool uniqOnly,std::vector<peak_content_t> peak_content, PeakIO *peaks, int min_length, std::string chrom, int smoothlen, std::vector<std::vector<double> > score_array_s)
{
/*"""Close the broad peak region, output peak boundaries, peak summit
and scores, then add the peak to peakIO object.

peak_content contains [start, end, treat_p, ctrl_p, list_scores]

peaks: a BroadPeakIO object

"""*/
    if(peak_content.size() ==0)
    {
        return false;
    }

    peak_content_t p_end = peak_content[peak_content.size() -1];
    peak_content_t p_start = peak_content[0];
    
    int peak_length = p_end.e - p_start.s ;
    if (peak_length < min_length) //if the peak is too small, reject it
    {
        return false;
    }
    
    std::vector<double> tlist_pileup;
    std::vector<double> tlist_control;
    std::vector<int> tlist_length;
    
    for(unsigned int i =0; i < peak_content.size(); i ++ )
    { // each position in broad peak
        peak_content_t p = peak_content[i];
    
        tlist_pileup.push_back( p.t_val );
        tlist_control.push_back( p.c_val );
        tlist_length.push_back( p.e - p.s );
    }
    //calculate p-value
   
    std::vector<double> tarray_pscore = __cal_pscore( tlist_pileup, tlist_control, uniqOnly );
    
   // std::vector<double> tarray_qscore = __cal_qscore( tlist_pileup, tlist_control );
   
    std::vector<double> tarray_fc     = __cal_FE    ( tlist_pileup, tlist_control );
   


    peaks->add( chrom,           // chromosome
               peak_content[0].s, // start
               peak_content[peak_content.size()-1].e, // end
               0, // summit position
               //mean_from_value_length( tarray_qscore, tlist_length ), // peak score
               0.0,
               mean_from_value_length( tlist_pileup, tlist_length ), // pileup
               mean_from_value_length( tarray_pscore, tlist_length ), // pvalue
               mean_from_value_length( tarray_fc, tlist_length ), // fold change
               //mean_from_value_length( tarray_qscore, tlist_length ) //qvalue
               0.0
               );
    
    return true;
}

