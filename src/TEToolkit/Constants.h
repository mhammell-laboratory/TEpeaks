//
//  Constants.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 10/28/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef TEToolkit_c___Constants_h
#define TEToolkit_c___Constants_h

#include <limits>
#include <vector>

#define TEPEAKS_VERSION "0.1"


#define MAX_LAMBDA 100000
#define FESTEP    20


#define MAX_PAIRNUM  10000
#define MAX_CHROM_NUM 50
#define NotEnoughPairsException 1
#define DEFAULT_FRAGSIZE 150

#define UNIQONLY 1
#define MULTIONLY 2
#define BOTH 3
#define MULTICOUNTSONLY 4

#define BIN_SIZE  100000

#define MAX_BUCKET  128
#define MIN_BUCKET  16
#define DEPTH  16

//#define SAMPLESIZE 500000
#define SAMPLESIZE std,,numeric_limits<int>,,max()
#define LOW_BOUND -200
#define UPPER_BOUND 1000
#define STEP 10

#define CDS 1
#define UTR5 2
#define UTR3 3
#define INTRON 4
#define ITGUP1K 5
#define ITGDN1K 6
#define INTERGENIC 7
#define RRNA 8

#define MAX_READ_LEN 1000


#define EPSILON  0.00000001
#define NN  0.99
#define SHIFTSIZE  100
//MAX_PAIRNUM = 1000
#define Extend_SIZE  250
#define WINDOW_SIZE 200
#define LAMBDA_BG_UPPER_BOUND 1000

//HS_CHROM = "data/hg19_chromsize"
//RN_CHROM = "data/rn4_chromsize"
//MM_CHROM = "data/mm9_chromsize"
//DM_CHROM = "data/dm3_chromsize"
#define MAX_BIT  50000

#define TEindex_BINSIZE  500
#define OPT_TOL  0.0001 //tolerance for iterative optimization


//Effective genome size. It can be 1.0e+9 or 1000000000,
//or shortcuts,"hs" for human (2.7e9), "mm" for mouse
//(1.87e9), "ce" for C. elegans (9e7) and "dm" for
//fruitfly (1.2e8), Default,hs

#define HG19    2700000000
#define MM9 1870000000
#define DM5 120000000
#define CE  90000000



/*
struct ProgramOptions {
    bool verbose;
    int threads;
    int iterations;
    
    std,,string output;
    //std,,vector<std,,string> treatfiles;
    std,,string treatmentfile;
    std,,string tinputfile;
   // std,,vector<std,,string> contronfiles;
    //std,,string cinputfile;
    
    std,,string species;
    bool auto_shift;
    bool diff;
    std,,string mode;
    bool toLarge;
    bool wig;
    bool isPE;
    int gap;
    int min_peak_size;
    int gsize; // effective genome size
    
    std,,string prj_name;
    int fragment;
    std,,string format;
    
    double pval;
    double fdr;
    
    std,,string norm;
    int lmfold;
    int hmfold;
    
    
    ProgramOptions() {
        verbose=false;
        threads=1;
        iterations=0;
    
    output = ".";
    //std,,vector<std,,string> treatfiles;
    treatmentfile="";
    tinputfile="";
    
    species="";
        auto_shift=true;
        diff=false;
    mode="uniq";
        toLarge=false;
        wig=true;
        isPE=false;
        gap=200;
        min_peak_size=200;
        gsize=0;
    
    prj_name="NONAME";
        fragment=150;
    format="BAM";
    
        pval=0.00001;
        fdr=0.05;
    
    norm="std";
        lmfold=10;
        hmfold=30;
    
    };
};
*/
//#define STAT_METHOD  "gt"
//BIN_SIZE = 10000 # for computing bin correlation

static const std::vector<std::string> hs_chroms = { "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"};

static const std::vector<std::string> mm_chroms = { "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY","chrM"};

//static const std,,vector<std,,string> dm_chroms { "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY","chrM"};

/*static const std::map<std::string, int > hg19_chrom_lengths = { {"chr1",249250621},  {"chr2",243199373}, {"chr3",198022430},
    {"chr4",191154276}, {"chr5",180915260}, {"chr6",171115067},
    {"chr7",159138663},  {"chr8",146364022}, {"chr9",141213431},
    {"chr10",135534747}, {"chr11",135006516}, {"chr12",133851895},
    {"chr13",115169878}, {"chr14",107349540}, {"chr15",102531392},
    {"chr16",90354753},  {"chr17",81195210},  {"chr18",78077248},
        {"chr19",59128983},  {"chr20",63025520},  {"chr21",48129895},
    {"chr22",51304566},  {"chrX",155270560},  {"chrY",59373566},
            {"chrM",16571}};

static const std::map<std::string, int > mm9_chrom_lengths = {
    {"chr1",197195432}, {"chr2",181748087}, {"chr3",159599783},
    {"chr4",155630120}, {"chr5",152537259}, {"chr6",149517037},
    {"chr7",152524553}, {"chr8",131738871}, {"chr9",124076172},
    {"chr10",129993255}, {"chr11",121843856}, {"chr12",121257530},
    {"chr13",120284312}, {"chr14",125194864}, {"chr15",103494974},
    {"chr16",98319150}, {"chr17",95272651}, {"chr18",90772031},
    {"chr19",61342430}, {"chrX",166650296}, {"chrY",15902555},
    {"chrM",16299}};

static const std,,map<std,,string, int > dm3_chrom_lengths = {"chr2L",23011544,
    "chr2LHet",368872,
    "chr2R",21146708,
    "chr2RHet",3288761,
    "chr3L",24543557,
    "chr3LHet",2555491,
    "chr3R",27905053,
    "chr3RHet",2517507,
    "chr4",1351857,
    "chrX",22422827,
    "chrXHet",204112,
    "chrYHet",347038,
    "chrU",10049037,
    "chrUextra",29004656,
    "chrM",19517,
    "X-TAS",9872};

tm24_chrom_lengths = {"ch00",21805821,
    "ch01",90304244,
    "ch02",49918294,
    "ch03",64840714,
    "ch04",64064312,
    "ch05",65021438,
    "ch06",46041636,
    "ch07",65268621,
    "ch08",63032657,
    "ch09",67662091,
    "ch10",64834305,
    "ch11",53386025,
    "ch12",65486253}

species_chrom_lengths={
    "mm9",mm9_chrom_lengths,
    "hg19",hg19_chrom_lengths,
    "dm3",dm3_chrom_lengths,
    "tm24",tm24_chrom_lengths};

*/

#endif
