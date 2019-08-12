//
//  Fwtrak.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 5/3/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____Fwtrak__
#define __TEToolkit_c____Fwtrak__

#include <stdio.h>
#include <vector>
#include <string>
#include <map>
//#include "Parser.h"

struct pos_t
{
    int start;
    int end;
};
class ShortRead
{
public :
    
    //unique reads only
    std::vector<std::map<int,std::vector<int> > > __locations;
    std::map<int,std::vector<pos_t> >  __locations_pe; //for PE: chromID -> poslist
    
    //multi reads only
    std::vector<std::map<int,std::vector<int> > > __locations_m;
    std::map<int,std::vector<pos_t> >  __locations_pe_m;
    
    
    bool __sorted;
    //std::vector<std::map<int,std::vector<int> > > __dup_locations;
   
    bool __dup_sorted;
    bool __destroyed;
    bool __dup_separated;
    double redundant_rate;
    //double sf; //scaling factor
    
    //std::vector<std::string> references;
    std::map<std::string,int> rlengths;
    
    //int total; //total number of tags
    int total_uniq; //peakmodel
    int total_multi;
    int total;
    int max_dup_tag;
    int dup_total;
    //public object annotation;
    //public object dups;
    //int fw;
    long length; // PE, sum of fragment length, unique reads only
    
    int tsize; //tag size
    int fraglength; //average frament length
    bool isPE;
    
    int average_template_length; //PE average fragment length
    //std::vector<std::string> chromNames;
    std::map<std::string,int> chromNameToId;
    std::map<int, std::string> chromIdToName;
    
    
    //ShortRead(std::string, opt_t );
    ShortRead();
    ~ShortRead();
    
    void clearMultiReads();
    //parser
    void inc_uniq_read();
    void inc_multi_read();
    void add_loc(std::string chrom, int fivepos, int strand,double w);
    void add_loc_pe(std::string chrom, int fivepos, int endpos,double w);
    void append_shortreads(ShortRead *other);
    int get_total();
    int get_total_uniq();
    void set_redundant_rate(double val);
    void clean_m();
    
    //tag/fragment length
    void set_taglength(int tsize);
    void set_fraglength(int fraglength);
    //reference chromosome lengths
    std::map<std::string,int> get_rlengths(); //parser
    void set_rlengths(std::map<std::string, int> reflengths); //parser
    
    //remove duplicates
    void sort();
    void separate_dups( int maxint = 1 );
    void separate_dups_pe (  int maxint  = 1);
    void separate_dups_multi(int maxint = 1); //EM first then remove dup

    //Peakmodel
    std::vector<int> get_locations_by_chr_uniq(int chrID,int strand);
    int get_num_of_chroms();
    std::vector<std::string> get_chrom_names();
    int numOfUniqTags_by_chr(int chrID,int strand);
    
   
    //void addback_dups();
    //void filter_dup ( int32_t maxnum = -1);
    //filter_dup_dryrun()
    
    
    //tuple extract_region_tags ( self, str chromosome, int32_t startpos, int32_t endpos );
    //compute_region_tags_from_peaks ( self, peaks, func, int window_size = 100, float cutoff = 5 );
    //refine_peak_from_tags_distribution ( self, peaks, int window_size = 100, float cutoff = 5 );
    
    //peak detect
    void sample_num (int samplesize, int seed = -1);
    void sample_percent (double percent, int seed = -1 );
    
    std::pair<std::vector<int>,std::vector<double> > pileup_a_chromosome ( std::string chrom, std::vector<int> ds, std::vector<double> scale_factor_s, bool uniqOnly, double baseline_value = 0.0, bool directional = true, int end_shift = 0 );
    //control
    std::pair<std::vector<int>,std::vector<double> > pileup_a_chromosome_c ( std::string chrom, std::vector<int> ds, std::vector<double> scale_factor_s, bool uniqOnly, double baseline_value = 0.0 );
    
    std::pair<std::vector<int>,std::vector<double> > pileup_a_chromosome_pe ( std::string chrom, std::vector<double> scale_factor_s, bool uniqOnly, double baseline_value = 0.0);
    
    
    void print_to_bed(std::string fname);
};


#endif /* defined(__TEToolkit_c____Fwtrak__) */
