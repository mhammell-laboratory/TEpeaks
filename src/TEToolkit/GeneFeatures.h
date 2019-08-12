//
//  GeneFeatures.h
//  TEToolkit_c
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c__GeneFeatures__
#define __TEToolkit_c__GeneFeatures__

#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>
#include <bitset>
#include <utility>

#include "IntervalTree.h"
//#include "Constants.h"
#include <limits>

typedef struct {
    std::string chrom;
    int start;
    int end;
} chr_ITV;

//extern "C" bool itv_comp(Interval, Interval);



typedef std::map<std::string, std::vector<std::string> > gene_exon_Dict ;
typedef gene_exon_Dict::iterator gene_exon_Dict_It;

typedef std::map<std::string,std::map<int, IntervalTree *> > chrom_itvTree_Dict ;
typedef chrom_itvTree_Dict::iterator chrom_itvTree_Dict_itr;
//typedef Dict::const_iterator gene_exon_Dict_It;


class Gene {
public:
    std::vector<std::pair<int,int> > exons;
    std::vector<std::pair<int,int> > cds;
    std::string id;
    std::string strand;
    int min_start ;
    int max_stop ;
    int gene_actual_len;
    int stop_codon_st;
    int stop_codon_end;
    int start_codon_st;
    int start_codon_end;
    std::vector<std::pair<int,int> > utr5;
    std::vector<std::pair<int,int> > utr3 ;
    std::vector<std::pair<int,int> > intron ;
    std::vector<std::pair<int,int> > itg1k ;
    Gene(std::string gid, std::string ss);
    //Gene& operator= (const Gene& gg);
    //Gene(const std::string gid, const std::string ss);
    ~Gene();
    void set_stop_codon(int st, int end);
    void add_cds(int st, int end);
    void add_exons(int st,int end);
    void get_others();
};

class GeneFeatures{

public:
    std::vector<std::string> features;
    std::vector<std::vector<int> > gene_percentile_list;
    
    std::vector<int> gene_starts;
    std::vector<int> gene_ends;
    std::vector<int> gene_lengths;

    int total_exon;
    
    chrom_itvTree_Dict cds_exon_idx_plus ;
    chrom_itvTree_Dict cds_exon_idx_minus ;
    
    GeneFeatures(std::string GTFfilename,std::string id_attribute);
    ~GeneFeatures();

    std::vector<std::string> getFeatures() ;
    std::string get_name(int g);
    
    int get_start(int g);
    int get_stop(int g);
    int get_numofgenes();
    int exist_in_percentile_list(int gene,int pos);

    std::map<int,int> Gene_annotation(std::string chrom, std::vector<std::pair<int,int> > itv_list, std::string strand,std::vector<int> * mapped_exons);
    
    std::vector<std::pair<int,int> > get_exons(std::string chrom,int read1_end,int read2_start,std::string strand);
    
private:
    
    void read_features(std::string gff_filename, std::string id_attribute) ;
    void build_tree(std::map<std::string, std::map<std::string,Gene> > temp_plus, std::map<std::string, std::map<std::string,Gene> > temp_minus);

};

/*extern "C" {
    void quick_sort(std::vector<Interval> &intervals, int first, int last);
    int pivot(std::vector<Interval> &intervals, int first, int last);
    int get_first(const std::pair<int, int>& p);
    
    int get_last(const std::pair<int, int>& p);
};*/

#endif /* defined(__BAMQC_0_5__GeneFeatures__) */
