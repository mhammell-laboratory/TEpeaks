//
//  PeakIO.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 6/10/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#include "PeakIO.h"
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <math.h>
#include "myLog.h"

using boost::algorithm::join;
using boost::lexical_cast;

char subpeak_letters(short i)
{
    if (i < 26)
    {
       
        return char(97 + i);
    }
    else
    {
        return subpeak_letters(i / 26) + char(97 + (i % 26));
    }
}

PeakIO::PeakIO(){
}
PeakIO::~PeakIO()
{
    
}


std::vector<PeakContent> PeakIO::get_data_from_chrom (std::string chrom)
{
    return peaks[chrom];
}

std::vector<std::string> PeakIO::get_chr_names ()
{
    std::vector<std::string> namelist;
    for(auto const& pairs: peaks)
        namelist.push_back(pairs.first);
    
    return namelist;
}

bool peakcontent_comp(PeakContent first, PeakContent second){
    return first.start < second.start ;
}

void PeakIO::sort ()
{
    //sort by position
    std::vector<std::string> chrs ;
    
    for (auto p : peaks){
        chrs.push_back(p.first);
    }
    std::sort(chrs.begin(),chrs.end());
    
    for ( auto chrom : chrs){
        std::sort(this->peaks[chrom].begin(),this->peaks[chrom].end(),peakcontent_comp);
    }
    
}


void PeakIO::add_PeakContent ( std::string chromosome, PeakContent pc )
{
    if (peaks.find(chromosome) != peaks.end())
    {
        peaks[chromosome].push_back(pc);
    }
    else {
        std::vector<PeakContent> plist;
        plist.push_back(pc);
        
        peaks.insert(std::pair<std::string,std::vector<PeakContent> >(chromosome,plist));
    }
        
}

void PeakIO::add (std::string chromosome, int start, int end, int summit ,
                  double peak_score, double pileup,
                  double pscore, double fold_change, double qscore)
{
    PeakContent pt ;
    
    pt.start = start;
    pt.end = end;
    pt.summit = summit;
    pt.score = peak_score;
    pt.pileup = pileup;
    pt.pscore = pscore;
    pt.qscore = qscore;
    pt.fc = fold_change;
    
    
    if (peaks.find(chromosome) == peaks.end())
    {
            std::vector<PeakContent> plist;
            plist.push_back(pt);
            
            peaks.insert(std::pair<std::string,std::vector<PeakContent> >(chromosome,plist));
    }
    else {
        peaks[chromosome].push_back(pt);
    }
    
}
void PeakIO::merge_peaks(PeakIO *pp)
{
    for(auto p : pp->peaks)
    {
        for (auto pt : p.second) {
            add(p.first,pt.start,pt.end,pt.summit,pt.score,pt.pileup,pt.pscore,pt.fc,pt.qscore);
        }
    }
}
//filter peaks by fold change
void PeakIO::filter_fc (double fc_low, double fc_up )
{
       // """Filter peaks in a given fc range.

       //If fc_low and fc_up is assigned, the peaks with fc in [fc_low,fc_up)
        
       // """
    
    std::map<std::string, std::vector<PeakContent> > new_peaks;
    
    if (fc_up != 0){
        for( auto p : peaks ){
            std::vector<PeakContent> plist;
            
            for (auto pt : p.second) {
                if (pt.fc >= fc_low and pt.fc < fc_up) {
                    plist.push_back(pt);
                }
            }
            if (plist.size() >0) {
                new_peaks.insert(std::pair<std::string, std::vector<PeakContent> > (p.first,plist) );
            }
            
        }
    }
    else{
        for( auto p : peaks ){
            std::vector<PeakContent> plist;
            
            for (auto pt : p.second) {
                if (pt.fc >=  fc_low) {
                    plist.push_back(pt);
                }
            }
            if (plist.size() > 0) {
                new_peaks.insert(std::pair<std::string,std::vector<PeakContent> >(p.first,plist) );
            }
            
        }
    }
    this->peaks = new_peaks;
    
}

//total number of peaks
int PeakIO::total ()
{
    int x = 0;
    for(auto p : peaks){
        x += p.second.size();
    }
    return x;
}

void PeakIO::write_to_xls (std::string fname, std::string narrowPeakFname, std::string name,std::string name_prefix, bool trackline,std::string score_column)
{

    std::ofstream fhd;
    
    fhd.open(fname,std::ofstream::out);
    
    if (peaks.size() > 0)
    {
        std::vector<std::string> strlist;
        strlist.push_back("chr");
        strlist.push_back("start");
        strlist.push_back("end");
        strlist.push_back("length");
        strlist.push_back("abs_summit");
        strlist.push_back("pileup");
        strlist.push_back("-log10(pvalue)");
        strlist.push_back("fold_enrichment");
        strlist.push_back("-log10(qvalue)");
        strlist.push_back("name");
        std::string joined = join(strlist,"\t");
        
        fhd <<  joined << std::endl;
    }
    else
    {
        fhd.close();
        return;
    }
    std::string peakprefix = name + "_" + name_prefix;

    std::ofstream fhd2;
    
    fhd2.open(narrowPeakFname,std::ofstream::out);

    if (trackline){
        //fprintf(fhd,"track type=narrowPeak name=\"%s\" description=\"%s\" nextItemButton=on\n", name, name);
        
        fhd2 << "track type=narrowPeak name=\"" << name << "\" description=\"" << name << "\" nextItemButton=on" << std::endl;
        
    }
    
    int n_peak = 0;
    int s = 0;
    
    for(auto p : peaks)
    {
        //debug("peaks.second size " + std::to_string(p.second.size()));
        std::map<int,std::vector<PeakContent> > ptlist_sortBy_end;
        
        for (auto pt : p.second) {
            if (ptlist_sortBy_end.find(pt.end) != ptlist_sortBy_end.end()) {
                ptlist_sortBy_end[pt.end].push_back(pt);
            }
            else {
                std::vector<PeakContent> ptmplist;
                ptmplist.push_back(pt);
                ptlist_sortBy_end.insert(std::pair<int,std::vector<PeakContent> >(pt.end,ptmplist));
            }
        }
        //debug("in write to xls ptlist_sortBy_end.size()" + std::to_string(ptlist_sortBy_end.size()));
        for (auto pp : ptlist_sortBy_end) {
            n_peak += 1;

            if (pp.second.size() >= 1){
                for( size_t i = 0; i < pp.second.size(); i ++ )
                {
                    PeakContent pt = pp.second[i];
                    std::string peakname = peakprefix + lexical_cast<std::string>(n_peak);
                    
                    fhd << p.first << "\t" << pt.start + 1 << "\t" << pt.end << "\t" << pt.end - pt.start + 1 << "\t" << pt.summit + 1 << "\t" << pt.pileup << "\t" << pt.pscore << "\t" << pt.fc << "\t" << pt.qscore << "\t" << peakname << std::endl;
                    
                    if (pt.summit == -1){
                        s = -1;
                    }
                    else{
                        s = pt.summit - pt.start;
                    }
                    double score = pt.score;
                    if (score_column == "pscore") {
                        score = pt.pscore;
                    }
                    if (score_column == "qscore") {
                        score = pt.qscore;
                    }
                    
                    
                 //   fprintf(fhd,"\t%d" , pt.summit +1); // summit position
                 //   fprintf(fhd,"\t%.2f", pt.pileup); // pileup height at summit
                 //   fprintf(fhd,"\t%.5f", pt.pscore); // -log10pvalue at summit
                 //   fprintf(fhd,"\t%.5f", pt.fc); // fold change at summit
                 //   fprintf(fhd,"\t%.5f", pt.qscore); // -log10qvalue at summit
                 //   fprintf(fhd,"\t%s", peakname);
                 //   fprintf("\n");
                    fhd2 << p.first << "\t" << pt.start << "\t" << pt.end << "\t" << peakname << "\t" << int(10*score) << "\t" << pt.fc << "\t" << pt.pscore << "\t" << pt.qscore << "\t" << s << std::endl;
                }
            }
        }
    }
    fhd.close();
    fhd2.close();
    return;
}

//these methods are very fast, specifying types is unnecessary
void PeakIO::write_to_bed (std::string fname, std::string name_prefix, std::string name,
                   std::string description , std::string score_column, bool trackline)
{
      /*  """Write peaks in BED5 format in a file handler. Score (5th
        column) is decided by score_column setting. Check the
        following list. Name column ( 4th column) is made by putting
        name_prefix together with an ascending number.

        Five columns are chromosome, peak start, peak end, peak name,
        and peak score.

        items in peak hash object:
      */
    //return _to_bed(name_prefix=name_prefix, name=name,
    //                        description=description, score_column=score_column,
    //                        print_func=fhd.write, trackline=trackline);
    std::string peakprefix = name + "_" + name_prefix;
    
    std::string desc = name + "_" + description;
    int n_peak = 0;
    
    std::string::size_type n = 0;
    while ( ( n = name.find( "\"", n ) ) != std::string::npos )
    {
        name.replace( n, 1, "\\\"" );
        n += 2;
    }
    
    n = 0;
    while ( ( n = desc.find( "\"", n ) ) != std::string::npos )
    {
        desc.replace( n, 1, "\\\"" );
        n += 2;
    }
    
    std::ofstream fhd;
    
    fhd.open(fname,std::ofstream::out);
    
    if (trackline){
        //fprintf(fhd,'track name="%s (peaks)" description="%s" visibility=1\n', name,desc);
        
        fhd << "track name=\"" << name << " (peaks)\" description=\"" << desc << "\" visibility=1\n" << std::endl;
        
    }
    for (auto p : peaks ){
        std::map<double,std::vector<PeakContent> > ptlist_sortBy_end;
        
        for (auto pt : p.second) {
            if (ptlist_sortBy_end.find(pt.end) != ptlist_sortBy_end.end()) {
                ptlist_sortBy_end[pt.end].push_back(pt);
            }
            else {
                std::vector<PeakContent> ptmplist;
                ptmplist.push_back(pt);
                ptlist_sortBy_end.insert(std::pair<double,std::vector<PeakContent> >(pt.end,ptmplist));
            }
        }
        for (auto pp : ptlist_sortBy_end) {
            n_peak += 1;
            
            if (pp.second.size() > 1) {
                for (size_t  i=0 ; i < pp.second.size(); i++ )
                {
                    PeakContent pt = pp.second[i];
                    double score = pt.score;
                    if (score_column == "pscore") {
                        score = pt.pscore;
                    }
                    if (score_column == "qscore") {
                        score = pt.qscore;
                    }
                    //fprintf(fhd,"%s\t%d\t%d\t%s%d%s\t%.5f\n", p.first,pt.start, pt.end,peakprefix,n_peak,subpeak_letters(i),score);
                    
                    fhd << p.first << "\t" << pt.start << "\t" << pt.end << "\t" << peakprefix << "\t" << n_peak << "\t" << subpeak_letters(i) << "\t" << score << std::endl;
                    
                }
            }
            else {
                PeakContent pt = pp.second[0];
                double score = pt.score;
                if (score_column == "pscore") {
                    score = pt.pscore;
                }
                if (score_column == "qscore") {
                    score = pt.qscore;
                }
                //fprintf(fdh,"%s\t%d\t%d\t%s%d\t%.5f\n", p.first,pt.start,pt.end,peakprefix,n_peak,score);
                
                fhd << p.first << "\t" << pt.start << "\t" << pt.end << "\t" << peakprefix << "\t" << n_peak << "\t" << score << std::endl;
            }
        }
    }
    
    fhd.close();
}

//these methods are very fast, specifying types is unnecessary

/*void PeakIO::append_candidate_to_bed (std::string chrom,PeakContent pt)
{
    
    std::ofstream ofs ;
    
    try {
        

        ofs.open (CandidatePeakFile, std::ofstream::out);
    
        if (ofs.is_open()) { //file already exists.
    
            ofs.close();
        
            ofs.open(CandidatePeakFile, std::ofstream::app);
            ofs << chrom << "\t" << pt.start << "\t" << pt.end << "\t" << "1" << "\t"  << pt.pileup << "\t" << pt.pscore << std::endl;
            ofs.close();
         }
        else { //first time write to the file.
            ofs << "chrom\tstar\tend\tpeakID\tpileup\tpscore\n";
            ofs << chrom << "\t" << pt.start << "\t" << pt.end << "\t" << "1" << "\t"  << pt.pileup << "\t" << pt.pscore << std::endl;
            ofs.close();
        }
    } catch (std::ofstream::failure e) {
        error("cannot output candidate peaks.");
        std::exit(1);
    }

}*/
void PeakIO::write_candidate_to_bed (std::string CandidatePeakFile)
{
    /*  """Write peaks in BED5 format in a file handler. Score (5th
     column) is decided by score_column setting. Check the
     following list. Name column ( 4th column) is made by putting
     name_prefix together with an ascending number.
     
     Five columns are chromosome, peak start, peak end, peak name,
     and peak score.
     
     items in peak hash object:
     */
    //return _to_bed(name_prefix=name_prefix, name=name,
    //                        description=description, score_column=score_column,
    //                        print_func=fhd.write, trackline=trackline);
    int n_peak = 0;
    std::ofstream ofs ;
    ofs.open (CandidatePeakFile, std::ofstream::out);
    info("candidate peak file name " + CandidatePeakFile);
    ofs << "chrom\tstar\tend\tpeakID\tpileup\tpscore\n";
    int i=0;
    info("candidate peaks size "  +  std::to_string(peaks.size()));
    for (auto p : peaks ){
        std::map<double,std::vector<PeakContent> > ptlist_sortBy_end;
        //debug("num of candidate peaks " + std::to_string(p.second.size()));
        
        for (auto pt : p.second) {
            //debug(std::to_string(i) + " pt.end " + std::to_string(pt.end));
            i += 1;
            if (ptlist_sortBy_end.find(pt.end) != ptlist_sortBy_end.end()) {
                ptlist_sortBy_end[pt.end].push_back(pt);
            }
            else {
                std::vector<PeakContent> ptmplist;
                ptmplist.push_back(pt);
                ptlist_sortBy_end.insert(std::pair<double,std::vector<PeakContent> >(pt.end,ptmplist));
            }
        }
        for (auto pp : ptlist_sortBy_end) {
            n_peak += 1;
            
            if (pp.second.size() > 1) {
                for (size_t i=0 ; i < pp.second.size(); i++ )
                {
                    PeakContent pt = pp.second[i];
                    
                    ofs << p.first << "\t" << pt.start << "\t" << pt.end << "\t" << n_peak << "\t"  << pt.pileup << "\t" << pt.pscore << std::endl;
                    //"%s\t%d\t%d\t%s%d%s\t%d\t%.5f\n", p.first,pt.start, pt.end,peakprefix,n_peak,subpeak_letters(i),pt.pileup,pt.pscore);
                }
            }
            else {
                PeakContent pt = pp.second[0];
                
                //fprintf(fhd,"%s\t%d\t%d\t%s%d\t%d\t%.5f\n", p.first,pt.start,pt.end,peakprefix,n_peak,pt.pileup,pt.pscore);
                
                ofs << p.first << "\t" << pt.start << "\t" << pt.end  << "\t" << n_peak << "\t" << pt.pileup << "\t" << pt.pscore << std::endl;
            }
        }
    }
    
    ofs.close();
}




