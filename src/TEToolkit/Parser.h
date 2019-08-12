//
//  Parser.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 5/3/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____Parser__
#define __TEToolkit_c____Parser__

//#include <boost/python.hpp>

#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include <pthread.h>
#include "htslib/sam.h"
#include <vector>
#include <map>

#include "Candidate_Peaks.h"
//#include "TEToolkit/Parse_BED.h"
#include "ShortRead.h"
#include "Constants.h"
#include "myLog.h"
//#include "Parser.h"

//using boost::python;

#if defined(MACOS)
//typedef pthread_mutex_t subread_lock_t;
#define pthread_spinlock_t pthread_mutex_t
#define pthread_spin_lock pthread_mutex_lock
#define pthread_spin_unlock pthread_mutex_unlock
#define pthread_spin_init(a, b) pthread_mutex_init(a, NULL)
#define pthread_spin_destroy(a) pthread_mutex_destroy(a)
#define strnlen(a,l) strlen(a)

#endif


struct opt_t
{
    int tsize;
    
    double gsize;
    double ratio; //ratio between IP and input treat/control
    double treat_scale;
    double ctrl_scale;
    int fragsize;
    
    //int extsize;
    int d;
    int scanwindow;
    int bandwidth;
    int shift;
    
    int lmfold;
    int umfold;
    
    int seed;
    int smalllocal;
    int largelocal;
    
    
    int threadNum;
    int numItr;
    int keepdupNum;
    
    bool log_pvalue;
    bool log_qvalue;
    double pval;
    double fdr;
    int pileup;
    double fe;
    double broadcutoff;
    double log_broadcutoff;
    
    std::string tfile;
    std::string cfile;
    std::string project_name;
    std::string data_outdir;
    std::string norm;
    std::string format;
    std::string modelR;
    std::string keepDuplicates;
    std::string species;
    
    bool PE_mode;
    bool onauto;
    //bool TEmode;
    bool tolarge;
    bool verbose;
    bool wig;
    bool broad;
    
    
    //bool nolambda;
    bool tocontrol;
    bool call_summits;
    bool downsample;

    std::string CandidatePeakfile;
    std::string peakxls;
    std::string peakNarrowPeak;
    bool trackline;
    std::string summitbed;
    std::vector<std::string> chromlist;
    
    opt_t() {
        verbose=false;
        threadNum=1;
        numItr=50;
        pileup = 20; //minium pileup required for peaks with multi-reads
        fe = 3; //minium fold enrichment for peaks with multi-reads
        
        data_outdir = ".";
        //std::vector<std::string> treatfiles;
        tfile="";
        cfile="";
        format ="BAM";
        bandwidth = 300;
        log_broadcutoff = 1;
        broadcutoff = 0.1;
        broad = false;
        
        species="hs";
        onauto=false;
        
        tolarge=false;
        wig=true;
        PE_mode=false;
        gsize= 0;
        keepdupNum=0;
        
        project_name="NONAME";
        keepDuplicates="1";
        fragsize=200;
        shift=0;
        
        pval=0.00001;
        fdr=0.05;
        ratio=1.0;
        treat_scale = 1.0;
        ctrl_scale = 1.0;
        
        norm="std";
        lmfold=5;
        umfold=50;
        smalllocal=1000;
        largelocal=10000;
        log_pvalue = false;
        log_qvalue = true;
        
    };

    
};

#define IS_PAIRED(bam) ((bam)->core.flag&BAM_FPAIRED)
#define IS_QCFAIL(bam) ((bam)->core.flag & BAM_FQCFAIL)

#define IS_PAIRED_AND_MAPPED(bam) (((bam)->core.flag&BAM_FPAIRED) && !((bam)->core.flag&BAM_FUNMAP) && !((bam)->core.flag&BAM_FMUNMAP))

#define IS_PROPERLYPAIRED(bam) (((bam)->core.flag&(BAM_FPAIRED|BAM_FPROPER_PAIR)) == (BAM_FPAIRED|BAM_FPROPER_PAIR) && !((bam)->core.flag&BAM_FUNMAP))

#define IS_UNMAPPED(bam) ((bam)->core.flag&BAM_FUNMAP)
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE)
#define IS_MATE_REVERSE(bam) ((bam)->core.flag&BAM_FMREVERSE)
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1)
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2)
#define IS_DUP(bam) ((bam)->core.flag&BAM_FDUP)




struct read_t
{
    std::string chrom;
    int start;
    int end;
    std::string strand;
    
};

struct peak_align_t
{
    std::vector<int> treat_multi_align_list;
    std::vector<int> ctrl_multi_align_list;
};

struct multi_read_t
{
    std::vector<int> pidlist;
    std::vector<int> startlist;
};

typedef std::map<int, multi_read_t> MULTI_READ_MAP;

template <class T>
struct read_pair_t {
    int readID;
    std::vector<T *> first_reads;
    std::vector<T *> second_reads;
} ;

template <class T>
struct thread_context_t
{
    unsigned short thread_id;
    pthread_t thread_object;
    pthread_spinlock_t cur_reads_lock;
    
    std::vector<read_pair_t<T> *> cur_reads;
    //std::vector<read_pair_t<T> *> * cur_reads;
    
    ShortRead * track;
    long fraglength;
    int n_frags;
    long taglength;
    int n_tags;
    
    std::vector<int> * peak_reads;
    MULTI_READ_MAP * multi_read_mapTo_pids;
    //MULTI_READ_MAP * multi_read_Align_start;
    std::map<int, std::vector<int> > * peak_uniqReads_startPos;
    
} ;

template <class T>
struct global_context_t
{
    unsigned short thread_number;
    int all_finished;
    int shiftsize;
    bool isPE;
    int flag;
    std::string format;
    std::vector<std::string> refnames;
    std::map<std::string,int> reflengths;
    std::vector<thread_context_t<T> *> thread_contexts;

    Candidate_Peaks * peakIdx;
} ;

template <class T>
struct arg_struct {
    global_context_t<T> * arg1;
    thread_context_t<T> * arg2;
};

extern "C" {
//void read_in_opts(std::string optfile,opt_t &options);

std::pair<ShortRead *, ShortRead *> read_aligmentFile(opt_t &options, int flag = UNIQONLY);

    void read_distribution(ShortRead * track, std::string inputFile,Candidate_Peaks * peakIdx,opt_t options,std::vector<double> * peak_reads_Prime, MULTI_READ_MAP * multi_read_mapTo_pids,std::map<int, std::vector<int> > * peak_uniqReads_startPos);
    
void read_distribution_BAM(ShortRead * track, std::string inputFile,Candidate_Peaks * peakIdx,opt_t options,std::vector<double> & peak_reads_Prime);
    
    void read_distribution_BED(ShortRead * track, std::string inputFile,Candidate_Peaks * peakIdx,opt_t options,std::vector<double> & peak_reads_Prime);
}

#endif /* defined(__TEToolkit_c____Parser__) */
