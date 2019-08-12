//
//  PeakIO.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 6/10/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____PeakIO__
#define __TEToolkit_c____PeakIO__

//#include <stdio.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include <fstream>


struct PeakContent
{
    int start;
    int end;
    int length;
    int summit;
    double height;
    double numTags;
    double score;
    double pileup;
    double pscore; //p-value
    double fc;
    double qscore;
    //std::string name;
};

class PeakIO
{
public:
    std::map<std::string,std::vector<PeakContent> > peaks;
    //std::string CandidatePeakFile;

    PeakIO();
    ~PeakIO();

    std::vector<PeakContent> get_data_from_chrom (std::string chrom);
    std::vector<std::string> get_chr_names ();
    void sort ();
    
    void filter_fc (double fc_low, double fc_up = 0 );
    int total();
    void add (std::string chromosome, int start, int end, int summit = 0,
               double peak_score=0, double pileup=0,
               double pscore=0, double fold_change=0, double qscore=0);
    void merge_peaks(PeakIO *p);
    void add_PeakContent ( std::string chromosome, PeakContent pc );
    //void set_CandidatePeakFile(std::string ofile);
    void write_to_bed (std::string fname, std::string name_prefix="peak_", std::string name="TEpeaks",
                       std::string description = "%s", std::string score_column="score", bool trackline=true);
    void write_to_xls (std::string fname, std::string narrowFile, std::string name="TEpeaks",std::string name_prefix="peak",bool trackline=true,std::string score_column="score");
    void write_candidate_to_bed (std::string CandidateFile);
    //void append_candidate_to_bed (std::string chrom,PeakContent pt);
    
//    void write_to_summit_bed (std::string fname, std::string name_prefix="peak_", std::string name="TEpeaks",
//                              std::string description="" , std::string score_column="score", bool trackline=true);
//    void write_to_narrowPeak (std::string fname,  std::string name="TEpeaks", std::string score_column="score", std::string name_prefix="peak", bool trackline=true);
};

#endif /* defined(__TEToolkit_c____PeakIO__) */
