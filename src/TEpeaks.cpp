//
//  TEpeaks.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 10/31/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
#include <thread>
#include <time.h>
#include <stdlib.h>
#include <math.h>


//#include "TEToolkit/Constants.h"
#include "Constants.h"
#include "TEToolkit/Parser.h"
//#include "narrow_TEpeaks.h"
using namespace std;
int run_narrow_TEpeaks_all(opt_t &opt);

bool ParseOptions_narrow(int argc, char **argv, opt_t& opt) {
    int verbose_flag = 0;
    int auto_flag = 0;
    int wig_flag = 0;
    int toLarge_flag = 0;
    int broad_flag = 0;
    //int diff_flag = 0;
    //int pe_flag = 0;
    
    //std::string norm = "sd";
    //double fdr = 0.05;
    
    
    
    const char *opt_string = "t:c:s:o:f:p:g:i:a:n:d:l:h:b:r:q:m:x:";
    
    static struct option long_options[] = {
        // long args
        {"verbose", no_argument, &verbose_flag, 1},
        {"auto", no_argument, &auto_flag, 1},
        {"wig", no_argument, &wig_flag, 1},
        {"broad",no_argument,&broad_flag, 1},
       // {"isPE",no_argument, &pe_flag,1},
        //{"diff", no_argument, &diff_flag, 1},
        //{"strand-specific", no_argument, &strand_flag, 1},
        {"toLarge", no_argument, &toLarge_flag, 1},
        {"treatment", required_argument, NULL, 't'},
        {"control", required_argument, NULL, 'c'},
        {"species", required_argument, NULL, 's'},
        {"outputdir", required_argument, NULL, 'o'},
        
        {"fraglen", required_argument, NULL, 'f'},
        {"pval",required_argument,NULL,'p'},
        {"gsize", required_argument, NULL, 'g'},
        {"prj_name", required_argument, NULL, 'n'},
        {"keepDup", required_argument, NULL,'d'},
        {"lmfold",required_argument,NULL,'l'},
        {"hmfold",required_argument,NULL,'h'},
        {"pileup", required_argument,NULL,'b'},
        {"iterations", required_argument, NULL, 'i'},

        {"fdr",required_argument,NULL,'q'},
        
        {"threads",required_argument, NULL, 'r'},

        {"shift", required_argument,NULL,'m'},
        
        {"fe", required_argument,NULL,'x'},
        {"broad-cutoff", required_argument, NULL, 'a'},
        {0,0,0,0}
    };
    
    int c;
    int option_index = 0;
    while (true) {
        c = getopt_long(argc,argv,opt_string, long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        
        switch (c) {
            case 0:
            
                break;
            case 'b' :{
                stringstream(optarg) >> opt.pileup;
                break;
            }
            case 'x' :{
                stringstream(optarg) >> opt.fe;
                break;
            }
           // case 'j' :{
           //     stringstream(optarg) >> opt.ratio;
           //     break;
           // }
        //    case 'k':{
         //       opt.format = optarg;
         //       break;
          //  }
            case 'm':{
                stringstream(optarg) >> opt.shift;
                break;
            }
            case 'c' : {
                stringstream(optarg) >> opt.cfile;
                break;
            }
            case 'g' : {
                stringstream(optarg) >> opt.gsize;
                break;
            }
            case 'a' : {
                stringstream(optarg) >> opt.broadcutoff ;
                opt.log_broadcutoff = log10(opt.broadcutoff)* (-1.0);
                break;
            }
            case 'r' : {
                stringstream(optarg) >> opt.threadNum;
                break;
            }
            case 'q' : {
                stringstream(optarg) >> opt.fdr;
                opt.log_qvalue = true;
                break;
            }
            case 'p' : {
                stringstream(optarg) >> opt.pval;
                opt.log_pvalue = true;
                
                break;
            }
            case 'l' : {
                stringstream(optarg) >> opt.lmfold;
                break;
            }
            case 'h' : {
                stringstream(optarg) >> opt.umfold;
                break;
            }
            case 'd' :{
                stringstream(optarg) >> opt.keepDuplicates;
                break;
            }
            case 't': {
                //IP sample files
                stringstream(optarg) >> opt.tfile;
                
                break;
            }
            case 'i': {
                
                stringstream(optarg) >> opt.numItr ;
                break;
            }
            case 's': {
                stringstream(optarg) >> opt.species;
                break;
            }
            case 'o': {
                //opt.data_outdir = optarg;
                stringstream(optarg) >> opt.data_outdir;
                break;
            }
            case 'n': {
                stringstream(optarg) >> opt.project_name;
                break;
            }
            case 'f': {
                
                stringstream(optarg) >> opt.fragsize;
                break;
            }
            
            default: break;
        }
    }
    
    
    if (verbose_flag) {
        opt.verbose = true;
    }
    
   if(broad_flag)
   {
       opt.broad = true;
   }
    if (auto_flag) {
        opt.onauto = true;
    }
    
    if (wig_flag) {
        opt.wig = true;
    }
    
 //   if (diff_flag) {
 //       opt.diff = true;
 //   }
    
    if (toLarge_flag) {
        opt.tolarge = true;
    }
    
    if (opt.tfile == "" || opt.cfile == "") {
        cout << "Please specify IP sample and input sample." << endl;
        exit(1);
    }
    if (opt.species != "hs" && opt.species != "mm" && opt.species!="ce" && opt.species != "dm" && opt.gsize <= 0 ) {
        
        cout << "Please specify species (hs,mm,ce,or dm are supported) or effective genome size using -g(--gsize) option!" << endl;
        exit(1);
    }
    
    if (opt.ratio <= 0) {
        opt.ratio = 1;
    }
   // if (opt.gsize >0) {
   //     opt.species = "";
   // }
    if (opt.species == "hs" ) {
        if (opt.gsize <=0) {
            opt.gsize = HG19;
        }
        
        opt.chromlist = hs_chroms;
        //cout << "chrom list size = " << hs_chroms.size() << endl;
    }
    if (opt.species == "mm") {
        if (opt.gsize <=0) {
            opt.gsize = MM9;
        }
        
        opt.chromlist = mm_chroms;
    }
    
    if(opt.species == "ce") {
        if (opt.gsize <=0) {
            opt.gsize = CE;
        }
    }
    if (opt.species == "dm") {
        if (opt.gsize <=0) {
            opt.gsize = DM5;
        }
    }
    
    if (opt.gsize < 10000) {
        cout << "Genome size is too small!" << endl;
        exit(1);
    }
    
    
    if (opt.keepDuplicates != "all" && opt.keepDuplicates !="auto") {
        int val = std::stoi(opt.keepDuplicates);
        
        opt.keepdupNum = val;
        if (val <= 0 || val >= 10) {
            cout << "Valid --keepDup option values are 'all', 'auto', all positive integer numbers!" << endl;
            exit(1);
        }
    }

/*    if (opt.norm != "std" && opt.norm != "bc") {
        cout << "Please specify normalization method supported by TEpeaks, 'std' or 'bc'!" << endl;
        exit(1);
    }*/
    
    //if user specified fragment length in a proper range, then use it directly, otherwise will run peak model to estimate fragment length.
    if (opt.fragsize < 50 || opt.fragsize > 500) {
        opt.onauto = true;
        opt.fragsize = 0;
        cout << "Fragment size is not properly set. TEpeaks will estimate fragment size automatically based on the data." << endl;
    }
    
    
    return true;
}

void PrintParam(opt_t opt){
    
    cout << "verbose = " << opt.verbose <<endl;
    cout << "num. of threads = " << opt.threadNum << endl;
    cout << "num. of iterations = " << opt.numItr << endl;
    
    cout << "output = " << opt.data_outdir << endl;
    //std::vector<std::string> treatfiles;
    /*for(size_t i =0 ; i < opt.treatfiles.size() ; i ++ )
    {
        cout << "treatment file = " << opt.treatfiles[i] << endl;
    }*/
    cout << "IP file = " << opt.tfile << endl;
    cout << "control file = " << opt.cfile << endl;
    //std::vector<std::string> contronfiles;
    //std::string cinputfile;
    
    cout << "species = " << opt.species << endl;
    cout << "gsize = " << opt.gsize << endl;
    if( opt.onauto)
    {
        cout << "auto = true" << endl;
    }
    else {
        cout << "auto = false" << endl;
    }
    
    if(opt.broad)
    {
        cout << "broad = true" << endl;
        cout << "broad-cutoff = " << opt.broadcutoff << endl;
    }
    
    //cout << "mode = " << opt.mode << endl;
    if(opt.tolarge)
    {
        cout << "tolarge = true" << endl;
    }
    else {
        cout << "tolarge = false" << endl;
    }
   // cout << "wig = " << opt.wig << endl;
    cout << "name = " << opt.project_name << endl;
    if(! opt.onauto)
    {
        cout << "fragment = " << opt.fragsize << endl;
    }
    cout << "bandwidth = " << opt.bandwidth << endl;
    
    cout << "pval = " << opt.pval << endl;
    cout << "fdr = " << opt.fdr << endl;
    
  //  cout << "norm = " << opt.norm << endl;
    cout << "lmfold = " << opt.lmfold << endl;
    cout << "hmfold = " << opt.umfold << endl;
    
    cout << "pileup = " << opt.pileup << endl;
    cout << "fe = " << opt.fe << endl;
    
    //cout << "is paired-end = " << opt.PE_mode<< endl;
}


void PrintVersion() {
    cout << "TEpeaks, version " << 	TEPEAKS_VERSION << endl;
}

void usage() {
    cout << "TEpeaks " << TEPEAKS_VERSION << endl << endl
    << "Usage: TEpeaks <CMD> [arguments] .." << endl << endl
    << "Where <CMD> can be one of:" << endl << endl
    << "    narrow         Call puntate peaks "<< endl
   // << "    broad        Call diffused peaks " << endl
   // << "    version       Prints version information"<< endl << endl
    << "Running TEpeaks <CMD> without arguments prints usage information for <CMD>"<< endl << endl;
}

void usage_narrow(bool valid_input = true) {
    if (valid_input) {
        
        cout << "TEpeaks  " << TEPEAKS_VERSION << endl
        << "Identifying transcription factor binding or histone modification sites" << endl << endl;
    }
    else {
        cout << "Invalid input!" << endl;
    }
    cout << "Usage: TEpeaks narrow [arguments] Alignment-files" << endl << endl
    << "Required arguments:" << endl
    << "-t, --treatment=STRING        IP sample ( BAM )" << endl
    << "-c, --control=files           Control IP samples ( BAM )"   << endl
    //<< "    --tinput=STRING          Input sample (BAM )" << endl
    //<< "    --cinput=STRING          control input sample" <<endl
    << "-o, --outputdir=STRING        Directory to write output to" << endl << endl

    
    << "Optional arguments:" << endl
  //  << "    --broad                  call broad peaks" << endl
  //  << "    --broad-cutoff           Cutoff for broad region. This option is not available unless --broad is set. If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff. (default: 0.1)" << endl
    //<< "    --auto                        Auto detect shiftsize for single-end reads (default: off)" << endl
    << "-f, --fraglen=INT             Fragment size (default: 200)" << endl
 //   << "    --ratio=DOUBLE            ratio between IP sample and input sample, IP/Input. By default, it's calculated from the data, but can also be set by user.  " << endl
    << "    --keepDup=STRING          How to deal with duplicate reads. The valid values are 'auto', 'all', or 1 (default: auto)" << endl
    << "    --shift=INT               Shift reads towards 3' end, if positive, or 5' end if negative. (default: 0)" << endl
    << "    --lmfold=INT              Lower limit of fold ratio against background to build model (default: 10)" << endl
    //<< "                              Model(default: 10)" << endl
    << "    --hmfold=INT              Higher limit of fold ratio against background to build model (default: 30)" << endl
    
    << "-n, --prjname=STRING          Name of the prject (default: NONAME)" << endl
    //<< "    --norm=STRING             normalization methods. sd (library size) or bc (bin" << endl
    //<< "                              correlation) (default: sd)" << endl
    << "-p, --pval=DOUBLE             P-value cutoff (default: 1e-5)" << endl
    << "    --fdr=DOUBLE              False discovery rate cutoff (default: 0.05)" << endl
   // << "-m, --mode=STRING             Running mode, e.g, multi or uniq (default: uniq)" << endl
    << "    --toLarge                 Scale library size to large sample (default: off)" << endl
    //<< "    --wig                     Output wiggle track (default: on)" << endl
    //<< "-d, --diff                    Differential peak analysis" << endl
    << "-s, --species=STRING          Species e.g., hs (Human hg19),  mm (Mouse mm9). (default: hs)" << endl
    << "-g, --gsize=INT               Effective genome size (default: human genome 2.7e9)" << endl
    
    << "    --threads=INT             Number of threads to use (default: 1)" << endl
    << "    --pileup=INT              The minuim pileup required for peaks with multi-reads (default: 20)" << endl
    << "    --fe=DOUBLE               The minuim fold enrichment required for peaks with multi-reads (default: 3)" << endl
    << "-i, --numItr=INT              Number of iterations (default: 50)" << endl;
    
}


//main function of TEpeaks. TEpeaks offers two subfunctions, narrow and broad to identify punctate peaks and diffused broad peaks respectively.
//plan to add differential peak analysis functions later.
int main(int argc, char *argv[]) {
    std::cout.sync_with_stdio(false);
    setvbuf(stdout, NULL, _IOFBF, 1048576);
    
    
    if (argc < 2) {
        usage();
        exit(1);
    } else {
        //clock_t tStart = clock();
        
        //a stucture to store all parameters
        opt_t opt;
        
        string cmd(argv[1]);
        
        if (cmd == "version") {
            PrintVersion();
            
        } else
        
        if (cmd == "narrow") {
            if (argc==2) {
                usage_narrow();
                return 0;
            }
            //parse paramters, return true if succeed.
            bool res = ParseOptions_narrow(argc-1,argv+1,opt);
            
            PrintParam(opt);
            
            if (!res) { //incorrect parameter settings
                cerr << endl;
                usage_narrow(false);
                
                exit(1);
            } else {
                
                // run the narrow Peak calling algorithm
                //take into account both unique reads and multi-reads.
                cout << "output = " << opt.data_outdir << endl;
                
                run_narrow_TEpeaks_all(opt);
                
                std::cout << "Done." << std::endl;
                
                cerr << endl;
            }
            //auto end_time(get_local_time());
            //clock_t tEnd = clock();
            //cerr << "Time spent: " << (double)(tEnd - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
            
        } else {
            cerr << "Error: invalid command " << cmd << endl;
            usage();
            exit(1);
        }
        
    }
    
    fflush(stdout);
    
    return 0;
}
