//
//  narrow_TEpeaks.cpp
//  TEToolkit_c++
//
//
//  Created by Ying Jin on 2/22/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#include "TEToolkit/EM_TEpeaks.h"

//#include <malloc.h>
#include <string>
#include <stdlib.h>
#include <pthread.h>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <ctime>
#include <sstream>
#include<dirent.h>

//#include <stdio>
#include <getopt.h>
#include <fstream>
#include <sys/stat.h>

//#include <boost/math/distributions/binomial.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/transformed.hpp>


//#include "IntervalTree.h"
//#include "TEToolkit/Candidate_Peaks.h"
//#include "TEToolkit/Constants.h"
#include "TEToolkit/zeroin.h"
#include "TEToolkit/myLog.h"
//#include "TEToolkit/Parser.h"
#include "TEToolkit/PeakDetect.h"
#include "TEToolkit/PeakModel.h"
//#include "TEToolkit/Output.h"
#include "TEToolkit/EM_TEpeaks.h"
#include "TEToolkit/ShortRead.h"

//using namespace boost::math::binomial;

using namespace boost::adaptors;
using namespace boost::algorithm;
using namespace boost;
//using namespace boost::lexical_cast;
using namespace std;

void model2r_script(PeakModel model, std::string filename, std::string name)
{
   // FILE *rfhd = fopen(filename,"w");
    std::ofstream rfhd ;
    
    rfhd.open (filename, std::ofstream::out);
    
    std::vector<float> norm_p ;
    std::vector<float> norm_m ;
    //#norm_s = [0]*w
    int sum_p = 0;
    int sum_m = 0;
    //#sum_s = sum(s)
    for (auto v : model.plus_line) {
        sum_p += v;
    }
    for (auto v : model.minus_line) {
        sum_m += v;
    }
    for (size_t i = 0 ; i < model.plus_line.size(); i++)
    {
        norm_p.push_back( 1.0 *(model.plus_line[i])*100/sum_p);
        norm_m.push_back( 1.0 *(model.minus_line[i])*100/sum_m);
    }
    
    rfhd << "# R script for Peak Model\n" << endl;
    
    std::string joined_norm_p = join( norm_p |transformed( static_cast<std::string(*)(float)>(std::to_string) ),
                                     ", " );
    std::string joined_norm_m = join( norm_m |transformed( static_cast<std::string(*)(float)>(std::to_string) ),
                                     ", " );
    std::string joined_ycorr = join( model.ycorr |transformed( static_cast<std::string(*)(float)>(std::to_string) ),
                                    ", " );
    std::string joined_xcorr = join( model.xcorr |transformed( static_cast<std::string(*)(float)>(std::to_string) ),
                                    ", " );
    std::string joined_alt_d = join( model.alternative_d |transformed( static_cast<std::string(*)(float)>(std::to_string) ),
                                    ", " );
    
    rfhd << "p <- c(" << joined_norm_p << ")" << endl;
    rfhd << "m <- c(" << joined_norm_m << ")" << endl;
    rfhd << "ycorr <- c(" << joined_ycorr << ")" << endl;
    rfhd << "xcorr <- c(" << joined_xcorr << ")" << endl;
    rfhd << "altd  <- c(" << joined_alt_d << ")" << endl;
    
    rfhd << "x <- seq.int((length(p)-1)/2*-1,(length(p)-1)/2)" << endl;
    
    rfhd << "pdf(' " << name << "_model.pdf',height=6,width=6)" << endl;
    
    rfhd << "plot(x,p,type='l',col=c('red'),main='Peak Model',xlab='Distance to the middle',ylab='Percentage')\n" << endl;
    rfhd << "lines(x,m,col=c('blue'))\n" << endl;
    rfhd << "legend('topleft',c('forward tags','reverse tags'),lty=c(1,1,1),col=c('red','blue'))\n" << endl;
    rfhd << "plot(xcorr,ycorr,type='l',col=c('black'),main='Cross-Correlation',xlab='Lag between + and - tags',ylab='Correlation')\n" << endl;
    rfhd << "abline(v=altd,lty=2,col=c('red'))\n" << endl;
    rfhd << "legend('topleft','alternative lag(s)',lty=2,col='red')\n" << endl;
    
    rfhd << "legend('right','alt lag(s) : " << joined_alt_d << "',bty='n')" << endl;
    
    rfhd << "dev.off()\n" << endl;
    
    
    rfhd.close();
}


int cal_max_dup_tags ( long genome_size, int tags_number, double p=1e-5 )
{
/*Calculate the maximum duplicated tag number based on genome
size, total unique tag number and a p-value based on binomial
distribution. Brute force algorithm to calculate reverse CDF no
more than MAX_LAMBDA(100000).*/
    
    /* inverts the binomial CDF. For lower tail only! */
    
    double b = 1.0 / genome_size ;
    
    return binomial_cdf_inv(1.0-p,tags_number,b);
    
}

//******consider both unique mappers and  multi-mappers***********
// part of the code is modified from MACS2.
int run_narrow_TEpeaks_all(opt_t &options)
{
    std::vector<double> IP_peakReads;
    std::vector<double> Input_peakReads;
    
    //a file to store candidate peak regions
    DIR * dp = opendir(options.data_outdir.c_str());
    
    if (dp!=NULL) {
        error("output dir already exists.");
        exit(1);
    }
    std::string cmd_creat_outputfolder = "mkdir " + options.data_outdir;
    
    
    info(cmd_creat_outputfolder);
    
    const int dir_err = system(cmd_creat_outputfolder.c_str());
    if (-1 == dir_err)
    {
        error("Error creating directory!n");
        exit(1);
    }
    
    options.CandidatePeakfile = options.data_outdir + "/" + options.project_name + "_candidatePeaks.txt";
    
    std::string tag = "tag" ;
    
    // STEP1 : read aligment files into ShortRead class
    info("#1 read alignment files ...");
    
    //read alignments of both unique reads and multi-reads
    //save unique reads into memory, for multi-reads， only count the numbers
    std::pair<ShortRead *, ShortRead *> treat_input_pair = read_aligmentFile(options,BOTH) ;
    //treatment sample total tag/fragment
    
    ShortRead *treat = treat_input_pair.first;
    ShortRead *control = treat_input_pair.second;
    
    if (!treat) {
        error("Reading in treatment file failed.");
        exit(1);
    }
    if (!control) {
        error("Reading in control file failed.");
        exit(1);
    }
    //paired-end sample don't need to buil peak model to estimate fragment length
    if(options.PE_mode)
    {
        tag = "fragment";
        options.onauto =false;
    }
    
    //tag size or fragment size (for paied-end sample) are estimated when parsing the alignment file.
    info("#1 " + tag + " size = " + std::to_string(options.tsize) );
    //get total mapped reads (both unique and multi-reads), multi-reads will be counted once only.
    int t0 = treat->get_total() ;
    int t0_u = treat->get_total_uniq();
    int t1 = t0;
    int t1_u = t0_u;
    int c0 = control->get_total(); //control sample
    int c0_u = control->get_total_uniq();
    int c1 = c0;
    int c1_u = c0_u;
    
    
    info("#1 total " + tag + "s in treatment: " + lexical_cast<std::string>(t0));
    int tm = treat->total_multi;
    info("#1 total multi " + tag + "s in treatment: " + lexical_cast<std::string>(tm));
    info("#1 total unique " + tag + "s in treatment: " + lexical_cast<std::string>(t0_u));
    int cm = control->total_multi;
    info("#1 total " + tag + "s in control: " + lexical_cast<std::string>(c0));
    info("#1 total multi " + tag + "s in control: " + lexical_cast<std::string>(cm));
    info("#1 total unique " + tag + "s in control: " + lexical_cast<std::string>(c0_u));
    
    int treat_max_dup_tags = 1;
    int control_max_dup_tags = 1;
    double treat_redundant_rate = 0;
    double control_redundant_rate = 0;
    
    // STEP2 : remove duplicates
    if (options.keepDuplicates != "all") {
        if (options.keepDuplicates == "auto") {
            //calculate max duplicate tags in single position based on binomial distribution
            info("#2 calculate max duplicate " + tag + "s in a single position based on binomial distribution ...");
            
            treat_max_dup_tags = cal_max_dup_tags(options.gsize,t0);
            //if ( treat_max_dup_tags == 0 ) { treat_max_dup_tags = 1;}
            info("#2 max duplicate in treatment based on binomial = " + lexical_cast<std::string>(treat_max_dup_tags));
            
            control_max_dup_tags = cal_max_dup_tags(options.gsize,c0);
            //if ( control_max_dup_tags == 0 ) { control_max_dup_tags = 1;}
            info("#2 max duplicate in control based on binomial = " + lexical_cast<std::string>(control_max_dup_tags));
        }
        else {
            //use user defined maximum tags
            info("#2 user defined maximum in treatment and control " + options.keepDuplicates);
            treat_max_dup_tags = lexical_cast<int>(options.keepDuplicates);
            control_max_dup_tags = lexical_cast<int>(options.keepDuplicates);
        }
        
        //remove duplicates
        info("#2 filter out redundant tags/fragments ...");
        if(options.PE_mode)
        {
            treat->separate_dups(treat_max_dup_tags);
            control->separate_dups(control_max_dup_tags);
        }
        else {
            treat->separate_dups(treat_max_dup_tags);
            control->separate_dups(control_max_dup_tags);
        }
        
        //using uniquely mapped reads to compute redundant rate
        t1_u = treat->get_total_uniq();
        c1_u = control->get_total_uniq();
        treat_redundant_rate = 1.0 * (t0_u - t1_u) / t0_u;
        control_redundant_rate =  1.0 * (c0_u - c1_u) / c0_u;
        
        treat->set_redundant_rate(treat_redundant_rate);
        control->set_redundant_rate(control_redundant_rate);
        
        //get total reads after removing duplicates
        t1 = treat->get_total();
        c1 = control->get_total();
        info("#2  total " + tag + "s after removing duplicates in treatment: " + lexical_cast<std::string>(t1));
        info("#2  total " + tag + "s after removing duplicates in control: " + lexical_cast<std::string>(c1));
    }
    
    //STEP3 : build peak model
    std::string joined_alt_d = "";
    
    //sort reads by position for each chromosome
    
    treat->sort();
    control->sort();
    
    
    if(!options.onauto) {
        info("#3 Skipped... building peak model to estimate fragment length.");
        
        if(options.PE_mode){
            options.d = options.tsize ;
        }
        else {
            options.d = options.fragsize ;
        }

        info("#3 Use " + lexical_cast<std::string>(options.d) + " as fragment length.");
        options.scanwindow = 2 * options.d  ;
    }
    else { //estimate fragment size based on data
        PeakModel pm(treat, options, MAX_PAIRNUM);
        info("#3 PeakModel finished!");
            // debug("#2  Summary Model:")
            // debug("#2   min_tags: %d" % (peakmodel.min_tags))
            // debug("#2   d: %d" % (peakmodel.d))
            // debug("#2   scan_window: %d" % (peakmodel.scan_window))
            
        info("#3 predicted fragment length is " + lexical_cast<std::string>(pm.d) + " bps");
            
        joined_alt_d = join( pm.alternative_d |
                                transformed( static_cast<std::string(*)(int)>(std::to_string) ),
                                ", " );
            
        info("#3 alternative fragment length(s) may be " + joined_alt_d);
            
        options.modelR = options.data_outdir + "/" +options.project_name + "_model.r";
        info("#3 Generate R script for peak model : " + options.modelR);
            
        model2r_script(pm,options.modelR,options.project_name);
        
        
        if (pm.NotEnoughPairs) {
            options.d = DEFAULT_FRAGSIZE;
            warn("#3 Not enough pairs found. Use default fragment length " +lexical_cast<std::string>(options.d) + " will be used as fragment length");
            
            options.scanwindow = 2 * options.d ;
            
        }
        else {
            options.d = pm.d;
            options.scanwindow = 2 * options.d;
            
            if (options.d <= 2 * options.tsize)
            {
                warn("#3 Since the d " + lexical_cast<std::string>(options.d) + " calculated from paired-peaks are smaller than 2*tag length, it may be influenced by unknown sequencing problem!");
                    
                warn("#3 Will use " + lexical_cast<std::string>(options.d) + " as fragment length. ");

                warn("#3 You may need to consider one of the other alternatives" );
                warn("#3 You can restart the process with -f --fragment XXX with your choice.");
                
            }
        }
        
    }
    
    //STEP4 : call peaks
    info("#4 Call peaks...");
//    if(options.nolambda){
//        info("#4 local lambda is disabled!");
//    }
    
    // decide tocontrol according to tolarge
    if(options.PE_mode){
        c1 = c1 * 2;
    }
    
      //  else {
            if (options.tolarge){
                if (t1 > c1)
                {
                    // treatment has more tags than control, since tolarge is
                    // true, we will scale control to treatment.
                    options.tocontrol = false;
                }
                else {
                    // treatment has less tags than control, since tolarge is
                    // true, we will scale treatment to control.
                    options.tocontrol = true;
                }
            }
            else {
                if (t1 > c1){
                    // treatment has more tags than control, since tolarge is
                    // false, we will scale treatment to control.
                    options.tocontrol = true;
                }
                else{
                    // treatment has less tags than control, since tolarge is
                    // false, we will scale control to treatment.
                    options.tocontrol = false;
                }
            }
     //   }

    if (options.shift > 0){
        info("#4 Sequencing ends will be shifted towards 3' by " + lexical_cast<std::string>(options.shift) + " bp(s) before calling peaks");
    }
    if (options.shift < 0) {
        info("#4 Sequencing ends will be shifted towards 5' by " + lexical_cast<std::string>(options.shift * -1) + " bp(s) before calling peaks" );
    }
    
    
    
    PeakDetect peakdetect( treat,control,&options);
    
    // call peaks on uniquely mapped reads only and output peaks
    peakdetect.call_peaks();
    
    info("#4 Output peaks...");
    std::string xlsfile = options.data_outdir + "/" + options.project_name + "_uniqpeaks.xls";
    
    
    std::string narrowFile = options.data_outdir + "/" + options.project_name +"_uniqpeaks.txt";
    //info("#4 Write peak in narrowPeak format file... " + options.peakNarrowPeak);
    //info("#4 Write peak in narrowPeak format file... " + narrowFile);
    //FILE *ofhd_bed = fopen( options.peakNarrowPeak, "w" );
    
    //peakdetect.peaks.write_to_narrowPeak (options.peakNarrowPeak, "peak_", options.project_name, score_column, options.trackline );
    std::string score_column = "pscore";
    //peakdetect.peaks.write_to_narrowPeak (narrowFile, "peak_", options.project_name, score_column, options.trackline );
    
    //peakdetect.peaks.write_to_xls(options.peakxls, options.project_name);
    peakdetect.peaks->write_to_xls(xlsfile,narrowFile,options.project_name);
    
    //peakdetect.peaks.write_to_narrowPeak(options);
    // call potential peak regions
    if(treat->total_multi > 0)
    {
    info("#4 Call potential peak regions for repetitive regions...");
   //peakdetect.peaks->set_CandidatePeakFile(options.CandidatePeakfile);
    peakdetect.candidate_peakregions();
    
    
    peakdetect.peaks->write_candidate_to_bed(options.CandidatePeakfile);
    //peakdetect.peaks->write_candidate_to_bed();
    
    treat->clean_m();
    control->clean_m();
    
    info("#4 Call peaks for repetitive regions...");
    run_EM_TEpeaks(options,treat,control,options.CandidatePeakfile);
           // filter out low FE peaks
    //peakdetect.peaks.filter_fc( fc_low = options.fecutoff );
    }
    
    delete treat;
    delete control;
    

    
    info("Done!");
    
    return 1;
}


