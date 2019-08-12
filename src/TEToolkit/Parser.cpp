//
//  narrow_TEpeaks.cpp
//  TEToolkit_c++
//
//
//  Created by Ying Jin on 2/22/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//


#include <string>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <ctime>
#include <sstream>

//#include <boost/iostreams/filter/zlib.hpp>

//#include <boost/iostreams/filter/gzip.hpp>
//#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/device/file.hpp>
#include <boost/algorithm/string/predicate.hpp>

//#include "IntervalTree.h"
//#include "TEToolkit/Candidate_Peaks.h"

//#include "TEToolkit/EMestimate_reads.h"
//#include "zeroin.h"
#include "myLog.h"
#include "Parser.h"
//#include "TEToolkit/PeakDetect.h"
//#include "TEToolkit/PeakModel.h"
//#include "ShortRead.h"

#include "htslib/sam.h"

unsigned int tick_time =3000;

/*bool isgz(std::string filename)
{
    //first, check that it exists
    std::ifstream in(filename);
    if(! in)
    {
        return false;
    }
    in.close();
    
    //
     If the file is not .gz and we try
     to read from it using boost's filtering_istream
     classes, an exception will be thrown.
     
     Therefore, catching the exception tells
     you the file is not .gz.
     //
    boost::iostreams::filtering_istream gzin;
    gzin.push(boost::iostreams::gzip_compressor());
    gzin.push(boost::iostreams::file_source(filename),std::ios_base::in | std::ios_base::binary);
    
    char c;
    try
    {
        gzin >> c;
    }
    catch ( boost::iostreams::gzip_error & e )
    {
        gzin.pop();
        return false;
    }
    return true;
}
*/

std::vector<std::pair<int, int> > fetch_exon(int st,uint32_t * cigar,uint32_t n_cigar,std::string format)
{
    //''' fetch exon regions defined by cigar. st must be zero based return list of tuple of (st, end)'''
    //match = re.compile(r'(\d+)(\D)')
    int chrom_st = st;
    if (format == "BAM") { chrom_st += 1;}
    
    std::vector<std::pair<int,int> >exon_bound;
    
    for (unsigned int i=0; i < n_cigar; i++){ //code and size
        if (bam_cigar_op(cigar[i])==BAM_CMATCH) { exon_bound.push_back(std::pair<int, int> (chrom_st,chrom_st + bam_cigar_oplen(cigar[i])-1));} //match
        if (bam_cigar_op(cigar[i])==BAM_CINS) {continue;} //insertion
        if (bam_cigar_op(cigar[i])==BAM_CDEL) { chrom_st += bam_cigar_oplen(cigar[i]);} //deletion
        if (bam_cigar_op(cigar[i])==BAM_CREF_SKIP) { chrom_st += bam_cigar_oplen(cigar[i]); } //gap
        if (bam_cigar_op(cigar[i])==BAM_CSOFT_CLIP) { chrom_st += bam_cigar_oplen(cigar[i]);}// soft clipping
    }
    return exon_bound;
}

int rightMost_pos(int st,uint32_t * cigar,uint32_t n_cigar,std::string format)
{
    //''' fetch exon regions defined by cigar. st must be zero based return list of tuple of (st, end)'''
    //match = re.compile(r'(\d+)(\D)')
    int chrom_st = st;
    if (format == "BAM") { chrom_st += 1;}
    
    //std::vector<std::pair<int,int> >exon_bound;
    
    for (unsigned int i=0; i < n_cigar; i++){ //code and size
        if (bam_cigar_op(cigar[i])==BAM_CMATCH) { chrom_st += bam_cigar_oplen(cigar[i])-1;} //match
        if (bam_cigar_op(cigar[i])==BAM_CINS) {continue;} //insertion
        if (bam_cigar_op(cigar[i])==BAM_CDEL) { chrom_st += bam_cigar_oplen(cigar[i]);} //deletion
        if (bam_cigar_op(cigar[i])==BAM_CREF_SKIP) { chrom_st += bam_cigar_oplen(cigar[i]); } //gap
        if (bam_cigar_op(cigar[i])==BAM_CSOFT_CLIP) { chrom_st += bam_cigar_oplen(cigar[i]);}// soft clipping
    }
    return chrom_st;
}

template <typename T>
void process_aligned_fragment(global_context_t<T> * global_context, thread_context_t<T> * thread_context,read_pair_t<T> * read_pair)
{
    bool multi_read = false;
    
    if (read_pair->first_reads.size() > 1 || read_pair->second_reads.size() >1) {
        multi_read = true;
    }
    
    //if (global_context->format=="BAM" || global_context->format == "SAM") {
        
        bam1_t * cur_read1 = NULL;
        bam1_t * cur_read2 = NULL;
        
        size_t alignement_cnt = 0;
        std::vector<std::pair<read_t, read_t> > multi_aligments;
        double w = 1.0;
        
        if (!multi_read ) { //unique read
            if (!global_context->isPE) { //SE
                cur_read1 = (bam1_t *)read_pair->first_reads[0];
                if((size_t)cur_read1->core.tid >= global_context->refnames.size()){
                    
                    error("Missing reference sequences with id " + std::to_string(cur_read1->core.tid));
                    std::exit(1);
                }
                else {
                    //reference sequence
                    std::string chrom = global_context->refnames[cur_read1->core.tid];
                    if (IS_REVERSE(cur_read1)) {
                        uint32_t *cigar = bam_get_cigar(cur_read1);
                        thread_context->track->add_loc(chrom,rightMost_pos(cur_read1->core.pos, cigar, cur_read1->core.n_cigar,global_context->format),1,w);
                        
                    }
                    else {
                        thread_context->track->add_loc(chrom,cur_read1->core.pos,0,w);
                    }
                    thread_context->taglength += cur_read1->core.l_qseq;
                    thread_context->n_tags += 1;
                    thread_context->track->inc_uniq_read();
                }
            }
            else {//PE
                
                if (read_pair->first_reads.size()==1) {
                    cur_read1 = (bam1_t *)read_pair->first_reads[0];
                    
                }
                if (read_pair->second_reads.size()==1) {
                    cur_read2 = (bam1_t *)read_pair->second_reads[0];
                }
                
                if(cur_read1 != NULL && !IS_UNMAPPED(cur_read1)&& !(cur_read1->core.flag & BAM_FQCFAIL) )
                {
                    std::string chrom1 = global_context->refnames[cur_read1->core.tid];
                    
                    //second read is not mapped, add first read only
                    if (cur_read2 == NULL || IS_UNMAPPED(cur_read2)||(cur_read2->core.flag & BAM_FQCFAIL)) {
                        uint32_t *cigar = bam_get_cigar(cur_read1);
                        
                        thread_context->track->add_loc_pe(chrom1,cur_read1->core.pos,rightMost_pos(cur_read1->core.pos, cigar, cur_read1->core.n_cigar,global_context->format),w);
                        
                    }
                    else {//both mates were mapped
                        std::string chrom2 = global_context->refnames[cur_read2->core.tid];
                        
                        int strand1 = IS_REVERSE(cur_read1)? -1 : 1;
                        int strand2 = IS_REVERSE(cur_read2)? -1 : 1;
                        
                        //same chromosome and different direction (considered proper pair)
                        if (chrom1 == chrom2 && (strand1 * strand2 ==-1) ) {
                            //tracking fragment length
                            thread_context->fraglength += std::abs(cur_read1->core.isize);
                            thread_context->n_frags += 1;
                            //using the left most position as start
                            int start = std::min(cur_read1->core.pos,cur_read2->core.pos);
                            //debug("add loc pe");
                            thread_context->track->add_loc_pe(chrom1,start,start + cur_read1->core.isize,w);
                            //debug("afater add loc pe");
                        }
                        else { //add bad mates separately
                            uint32_t *cigar = bam_get_cigar(cur_read1);
                            //debug("add loc pe bad mate1");
                            thread_context->track->add_loc_pe(chrom1,cur_read1->core.pos,rightMost_pos(cur_read1->core.pos, cigar, cur_read1->core.n_cigar,global_context->format),w*0.5);
                            //debug("afater add loc pe bad mate1");
                            cigar = bam_get_cigar(cur_read2);
                            //debug("add loc pe bad mate2");
                            thread_context->track->add_loc_pe(chrom2,cur_read2->core.pos,rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format),w*0.5);
                            //debug("afater add loc pe bad mate2");
                        }
                    }//both ends mapped
                    thread_context->track->inc_uniq_read();
                }//good alignment
                
                else { //read1 is not mapped
                    if(cur_read2 != NULL && !IS_UNMAPPED(cur_read2) && !(cur_read2->core.flag & BAM_FQCFAIL) )
                    {
                        std::string chrom = global_context->refnames[cur_read2->core.tid];
                        uint32_t *cigar = bam_get_cigar(cur_read2);
                    
                        thread_context->track->add_loc_pe(chrom,cur_read2->core.pos,rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format),w);
                        thread_context->track->inc_uniq_read();
                    }
                    
                }//end read1 is not mapped
                
            }//end PE
            
        }//end uniq-read
        else {
            //========================multi-read=========================
        //loop over aligments of the first segment
        for (size_t k =0;  k < read_pair->first_reads.size() ; k++) {
            
            alignement_cnt += 1;
            
            read_t algn1;
            read_t algn2;
            
            algn1.chrom = "";
            algn2.chrom = "";
            algn1.start = -1;
            algn2.start = -1;
            algn1.end = -1;
            algn2.end = -1;
            algn1.strand = ".";
            algn2.strand = ".";
            
            //get the first alignment of the first segment
            cur_read1 = (bam1_t *)read_pair->first_reads[k];
            
            if (read_pair->second_reads.size() > 0) {
                //the first alignment of the second segment
                cur_read2 = (bam1_t *)read_pair->second_reads[k];
            }
            
            if(cur_read1 != NULL && !IS_UNMAPPED(cur_read1) && !(cur_read1->core.flag & BAM_FQCFAIL)){
                    //chomosome appear in the header
                    if((size_t)cur_read1->core.tid >= global_context->refnames.size()){
                    
                        error("Missing reference sequences with id " + std::to_string(cur_read1->core.tid));
                        std::exit(1);
                    }
                    else {
                        //reference sequence
                        algn1.chrom = global_context->refnames[cur_read1->core.tid];
                    }
                    algn1.strand = IS_REVERSE(cur_read1) ? "-" : "+";
                    //left most position
                    algn1.start = cur_read1->core.pos;
                    
                    uint32_t *cigar = bam_get_cigar(cur_read1);
                    //save the right most position
                    algn1.end = rightMost_pos(cur_read1->core.pos, cigar, cur_read1->core.n_cigar,global_context->format);
                    
                    bam_destroy1(cur_read1);
            }
            
            if(cur_read2 != NULL && ! IS_UNMAPPED(cur_read2)&& !(cur_read2->core.flag & BAM_FQCFAIL)){
                    //uint32_t *cigar = bam_get_cigar(cur_read2);
                    
                    algn2.strand =  IS_REVERSE(cur_read2) ? "-" : "+";
                    
                    if((size_t)cur_read2->core.tid >= global_context->refnames.size()){
                        
                        error("Missing reference sequences with id " + std::to_string(cur_read2->core.tid));
                        std::exit(1);
                    }
                    else {
                        algn2.chrom = global_context->refnames[cur_read2->core.tid];
                    }

                    algn2.start = cur_read2->core.pos;
                    uint32_t *cigar = bam_get_cigar(cur_read2);
                    //save the right most position
                    algn2.end = rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format);
                    
                    bam_destroy1(cur_read2);
                
            }
            
            if (algn1.chrom != "" || algn2.chrom != "") {
                multi_aligments.push_back(std::pair<read_t, read_t> (algn1,algn2));
            }

        }//end loop cur_read1
        
        //there are remaining alignments of the second segment, e.g., only the second segment is mapped
        if (alignement_cnt < read_pair->second_reads.size()) {
            
            read_t algn2;
            algn2.chrom = "";
            algn2.start = -1;
            algn2.end = -1;
            algn2.strand = ".";
            
            for (size_t i = alignement_cnt; i < read_pair->second_reads.size(); i++) {
                
                read_t algn2;
                read_t algn1;
                
                algn2.chrom = "";
                algn2.start = -1;
                algn2.end = -1;
                algn2.strand = ".";
                
                algn1.chrom = "";
                algn1.start = -1;
                algn1.end = -1;
                algn1.strand = ".";
                
                cur_read2 = (bam1_t *)read_pair->second_reads[i];
                if(cur_read2 != NULL && ! IS_UNMAPPED(cur_read2)&& !(cur_read2->core.flag & BAM_FQCFAIL)){
                    if((size_t)cur_read2->core.tid >= global_context->refnames.size()){
                        error("Missing reference sequences with id " + std::to_string(cur_read2->core.tid));
                        std::exit(1);
                    }
                    algn2.chrom = global_context->refnames[cur_read2->core.tid];
                    algn2.start = cur_read2->core.pos;
                    uint32_t *cigar = bam_get_cigar(cur_read2);
                    //save the right most position
                    algn2.end = rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format);
                
                    algn2.strand =  IS_REVERSE(cur_read2) ? "-" : "+";
                
                    bam_destroy1(cur_read2);
                
                    if (algn1.chrom != "" || algn2.chrom != "") {
                        multi_aligments.push_back(std::pair<read_t, read_t> (algn1,algn2));
                    }
                }
            }//end loop cur_read2
            
        }// end if more cur_read2
            
        if (multi_aligments.size() >0) {
            thread_context->track->inc_multi_read();
        }
        if (global_context->flag !=MULTICOUNTSONLY ) {
                
        //weight 1/n
        w = 1.0 / multi_aligments.size();
            
        for (size_t k =0 ; k < multi_aligments.size(); k++) {
                
                read_t r1 = multi_aligments[k].first;
                read_t r2 = multi_aligments[k].second;
                
                if (r1.chrom != "" && r2.chrom !="") {
                    //proper pairs
                    if (r1.chrom == r2.chrom && (r1.strand != r2.strand)) {
                        int s = std::min(r1.start,r2.start);
                        int e = std::max(r1.end,r2.end);
                        
                        thread_context->track->add_loc_pe(r1.chrom,s,e,w);
                    }
                    /*else { //bad pairs add individual segments
                        thread_context->track->add_loc_pe(r1.chrom,r1.start,r1.end,w*0.5);
                        thread_context->track->add_loc_pe(r2.chrom,r2.start,r2.end,w*0.5);
                    }*/
                    
                }//both mate mapped
                else{
                    if (global_context->isPE) { //PE
                        //debug("before add loc " + std::to_string(w));
                        if (r1.chrom != "") {
                            thread_context->track->add_loc_pe(r1.chrom,r1.start,r1.end,w);
                        }
                        if (r2.chrom != "") {
                            thread_context->track->add_loc_pe(r2.chrom,r2.start,r2.end,w);
                        }
                    }
                    else {//SE
                        //debug("before add loc " + std::to_string(w));
                        if (r1.strand == "-") {
                            thread_context->track->add_loc(r1.chrom,r1.end,1,w);
                        }
                        else {
                            thread_context->track->add_loc(r1.chrom,r1.start,0,w);
                        }
                    }//end SE
                }//one mate mapped
                
            }//end loop
            
        }//end if
        
        }//end multi-read
        delete read_pair;
    
    
}

//parse multi-reads only
template <typename T>
void process_aligned_fragment_multi(global_context_t<T> * global_context,thread_context_t<T> * thread_context,read_pair_t<T> * read_pair)
{
    //debug("in process_aligned_fragment_multi");
    //if (global_context->format=="BAM" || global_context->format == "SAM") {
        
        bam1_t * cur_read1 = NULL;
        bam1_t * cur_read2 = NULL;
        int readID = read_pair->readID;
        std::string strand1 = ".";
        std::string strand2 = ".";
        std::vector<std::pair<int, int> > exons ;
        //std::string qname = "";
        
        std::string chrom = "";
        
        size_t alignement_cnt = 0;
        bool multi_read = false;
        
        if(read_pair->first_reads.size() >1 || read_pair->second_reads.size() > 1)
        {
            multi_read = true;
            //debug("found a multi-read " + std::to_string(readID));
        }
        
        //loop through each alignment
        for (size_t k =0;  k < read_pair->first_reads.size() ; k++) {
            int chrom1_id = -1;
            int chrom2_id = -1;
            alignement_cnt += 1;
            
            cur_read1 = (bam1_t *)read_pair->first_reads[k];
            
            if (read_pair->second_reads.size() > 0) {
                cur_read2 = (bam1_t *)read_pair->second_reads[k];
            }
            
            int pid1 = -1;
            int pid2 = -1;
            int pid = -1;
            int right_most_pos1 = -1;
            
            int right_most_pos2 = -1;
            int start_pos = -1;
            
            if(cur_read1 != NULL && !IS_UNMAPPED(cur_read1) && !(cur_read1->core.flag & BAM_FQCFAIL)){
                    uint32_t *cigar = bam_get_cigar(cur_read1);
                    //char * name = bam_get_qname(cur_read1);
                    chrom1_id  =  cur_read1->core.tid;
                    
                    //qname = std::string(name);
                    
                    if((size_t)cur_read1->core.tid < global_context->refnames.size()){
                        chrom = global_context->refnames[cur_read1->core.tid];
                    }
                    else {
                        std::cout << (size_t)cur_read1->core.tid  << "\t" << "missing reference sequences." << std::endl;
                        std::exit(1);
                    }
                
                    strand1 = IS_REVERSE(cur_read1) ? "-" : "+";
                    
                    exons = fetch_exon(cur_read1->core.pos, cigar,cur_read1->core.n_cigar, global_context->format);
                    
                    pid1 = global_context->peakIdx->get_ovp_peaks(chrom,exons,global_context->shiftsize,strand1);
                    if(strand1 == "+")
                    {
                        start_pos = cur_read1->core.pos;
                    }
                    else {
                        right_most_pos1 = rightMost_pos(cur_read1->core.pos, cigar, cur_read1->core.n_cigar,global_context->format);
                    }
                    
                    bam_destroy1(cur_read1);
                    
            }// end if(cur_read1 != NULL )
            
            if(cur_read2 != NULL && ! IS_UNMAPPED(cur_read2)&& !(cur_read2->core.flag & BAM_FQCFAIL)){
                    uint32_t *cigar = bam_get_cigar(cur_read2);
                    
                    strand2 =  IS_REVERSE(cur_read2) ? "-" : "+";
                    
                    chrom2_id = cur_read2->core.tid;
                    
                    if((size_t)cur_read2->core.tid < global_context->refnames.size()){
                        chrom = global_context->refnames[cur_read2->core.tid];
                    }
                    else {
                        std::cout << (size_t)cur_read2->core.tid << "\t" <<"missing reference sequences." << std::endl;
                        std::exit(1);
                    }
                    exons = fetch_exon(cur_read2->core.pos, cigar,cur_read2->core.n_cigar, global_context->format);
                    
                    pid2 = global_context->peakIdx->get_ovp_peaks(chrom,exons,global_context->shiftsize,strand2);
                    
                    if(strand2 == "-")
                    {
                        right_most_pos2 =rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format);
                    }
                    else {
                        
                        start_pos = cur_read2->core.pos;
                    }
                    bam_destroy1(cur_read2);
                
            } //end if(cur_read2 != NULL )
            
            //PE reads and both ends mapped
            if (chrom1_id != -1 && chrom2_id != -1){
                //discordantly mapped pairs
                if (chrom1_id != chrom2_id ){
                    continue;
                }
                if( pid1 != pid2) //overlap with different peaks, need to be resolved in the next version
                {
                    continue;
                }
                
            }
            //resolve overlapped candidate peak id
            if(pid1 == -1 && pid2  != -1)
            {
                pid = pid2;
                if(start_pos == -1)
                {
                    start_pos = right_most_pos2 - 2 * global_context->shiftsize;
                    
                }
            }
            if(pid2 == -1 && pid1 != -1)
            {
                pid = pid1;
                if(start_pos == -1)
                {
                    start_pos = right_most_pos1 - 2 * global_context->shiftsize;
                }
            }
            if(start_pos < 0)
            {
                start_pos = 0;
            }
            
            //overlap with a candidate peak
            if (pid != -1) {
                if (multi_read) {
                    //save peak ids mapped by the multi-read, and the start position
                    MULTI_READ_MAP::iterator itr;
                    itr = thread_context->multi_read_mapTo_pids->find(readID);
                    
                    if (itr != thread_context->multi_read_mapTo_pids->end()) {
                        (* thread_context->multi_read_mapTo_pids)[readID].pidlist.push_back(pid);
                        //keep multi_read starting position
                        (* thread_context->multi_read_mapTo_pids)[readID].startlist.push_back(start_pos);
                        
                    }
                    else {
                        multi_read_t m ;
                        
                        m.pidlist.push_back(pid); // = pidlist;
                        m.startlist.push_back(start_pos);
                        thread_context->multi_read_mapTo_pids->insert(std::pair<int,multi_read_t > (readID,m));
                       
                    }
                }//end multi-read
                
            else { //uniq read, save start position for the corresponding peak id
                    thread_context->peak_reads->at(pid) += 1.0;
                
                    std::map<int,std::vector<int> >::iterator itr;
                    itr = thread_context->peak_uniqReads_startPos->find(pid);
                    
                    if (itr != thread_context->peak_uniqReads_startPos->end()) {
                       (* thread_context->peak_uniqReads_startPos)[pid].push_back(start_pos);
                    }
                    else {
                        
                        std::vector<int> startlist;
                        startlist.push_back(start_pos);
                        thread_context->peak_uniqReads_startPos->insert(std::pair<int, std::vector<int> >(pid,startlist));
                    }
            } //end uniq read
            }//end pid
        }
        
        if (alignement_cnt < read_pair->second_reads.size()) {
            for (size_t i = alignement_cnt; i < read_pair->second_reads.size(); i++) {
                cur_read2 = (bam1_t *)read_pair->second_reads[i];
                int pid = -1;
                
                int right_most_pos2 = -1;
                int start_pos = -1;
                
                //int chrom2_id = -1;
                if(! IS_UNMAPPED(cur_read2)&& !(cur_read2->core.flag & BAM_FQCFAIL)){
                    uint32_t *cigar = bam_get_cigar(cur_read2);
                    
                    strand2 =  IS_REVERSE(cur_read2) ? "-" : "+";
                    
                    //chrom2_id = cur_read2->core.tid;
                    
                    if((size_t)cur_read2->core.tid < global_context->refnames.size()){
                        chrom = global_context->refnames[cur_read2->core.tid];
                    }
                    else {
                        std::cout << (size_t)cur_read2->core.tid << "\t" <<"missing reference sequences." << std::endl;
                        std::exit(1);
                    }
                    exons = fetch_exon(cur_read2->core.pos, cigar,cur_read2->core.n_cigar, global_context->format);
                    
                    pid = global_context->peakIdx->get_ovp_peaks(chrom,exons,global_context->shiftsize,strand2);
                    
                    if(strand2 == "-")
                    {
                        right_most_pos2 =rightMost_pos(cur_read2->core.pos, cigar, cur_read2->core.n_cigar,global_context->format);
                        start_pos = right_most_pos2 - 2 * global_context->shiftsize;
                        if(start_pos <0)
                        {
                            start_pos = 0 ;
                        }
                    }
                    else {
                        
                        start_pos = cur_read2->core.pos;
                    }
                    
                    bam_destroy1(cur_read2);
                }//cur_read2 is not null
                
                if (pid != -1) {
                    if (multi_read) {
                        //std::map<int,std::vector<int> >::iterator itr;
                        MULTI_READ_MAP::iterator itr;
                        itr = thread_context->multi_read_mapTo_pids->find(readID);
                        
                        if (itr != thread_context->multi_read_mapTo_pids->end()) {
                            (* thread_context->multi_read_mapTo_pids)[readID].pidlist.push_back(pid);
                            (* thread_context->multi_read_mapTo_pids)[readID].startlist.push_back(pid);
                            //keep multi_read starting position
                            //(* thread_context->multi_read_Align_start)[readID].push_back(start_pos);
                        }
                        else {
                            multi_read_t m;
                            //std::vector<int> pidlist ;
                            m.pidlist.push_back(pid);
                            m.startlist.push_back(start_pos);
                            thread_context->multi_read_mapTo_pids->insert(std::pair<int,multi_read_t > (readID,m));
                           
                        }
                        
                        
                        
                    }
                    else { //uniq read
                        
                        //    (*thread_context->peak_reads)[pid] += 1;
                        thread_context->peak_reads->at(pid) += 1.0;
                        
                        std::map<int,std::vector<int> >::iterator itr;
                        itr = thread_context->peak_uniqReads_startPos->find(pid);
                        
                        if (itr != thread_context->peak_uniqReads_startPos->end()) {
                            (* thread_context->peak_uniqReads_startPos)[pid].push_back(start_pos);
                        }
                        else {
                            
                            std::vector<int> startlist;
                            startlist.push_back(start_pos);
                            thread_context->peak_uniqReads_startPos->insert(std::pair<int, std::vector<int> >(pid,startlist));
                        }
                      
                    }//end uniq read
                    
                }//end pid
                
            }//loop all remaining cur_read2
        }//end if more cur_read2 alignments
        
        delete read_pair;

}


template <typename T>
void* worker(void * vargs)
{
    struct arg_struct<T> * args = (struct arg_struct<T> *) vargs;
    thread_context_t<T> * thread_context = args->arg2;
    global_context_t<T> * global_context = args->arg1;
    delete args;
    
    while (1){
        //std::pair<bam1_t *, bam1_t *> cur_read_pair;
        read_pair_t<T>  * cur_read_pair;
        while(1){
            int is_retrieved = 0;
            pthread_spin_lock(&thread_context->cur_reads_lock);
            
            if(thread_context->cur_reads.size()>0){
                cur_read_pair = thread_context->cur_reads.back();
                thread_context->cur_reads.pop_back();
                
                is_retrieved = 1;
            }
            pthread_spin_unlock(&thread_context->cur_reads_lock);
            if(global_context->all_finished && !is_retrieved) return NULL;
            
            if(is_retrieved) break;
            else usleep(tick_time);
        }
        if (global_context->flag == MULTIONLY) {
            process_aligned_fragment_multi<T>(global_context,thread_context,cur_read_pair);
        }
        else {
            process_aligned_fragment<T>(global_context,thread_context,cur_read_pair);
        }
        
        
    }
}


//start threads
template <typename T>
void init_thread(global_context_t<T> * global_context,unsigned short threadNumber)
{
    global_context->thread_number = threadNumber;
    
    int numOfpeaks = 0;
    if (global_context->peakIdx != nullptr) {
        numOfpeaks = global_context->peakIdx->get_numofpeaks();
    }
   
    
    if(threadNumber >1){
        for(int i=0;i<threadNumber;i++)
        {
            thread_context_t<T> *th_contx = new thread_context_t<T>();
            pthread_spin_init(&th_contx->cur_reads_lock, PTHREAD_PROCESS_PRIVATE);
            th_contx->thread_id = i;
            // th_contx->cur_reads = new std::vector<std::pair<bam1_t *,bam1_t *> >() ;
            
           //// th_contx->cur_reads = new std::vector<read_pair_t<T> *> ();
            th_contx->track = new ShortRead();
            
            if (numOfpeaks > 0) {
                //debug("init threads init peak_reads");
                th_contx->peak_reads = new std::vector<int>();
            
                for (int j=0; j < numOfpeaks; j++) {
                    th_contx->peak_reads->push_back(0);
                }
                th_contx->multi_read_mapTo_pids = new MULTI_READ_MAP ();
                th_contx->peak_uniqReads_startPos = new std::map<int, std::vector<int> > ();
            }
        
            global_context->thread_contexts.push_back(th_contx);
        }
        
        for(int i=0;i<threadNumber;i++)
        {
            struct arg_struct<T> * args = new arg_struct<T>();
            args->arg1 = global_context;
            args->arg2 = global_context->thread_contexts[i];
            //int ret = pthread_create(&global_context -> thread_contexts[i].thread_object, NULL, worker, (void *)&args);
            pthread_create(&global_context -> thread_contexts[i]->thread_object, NULL, worker<T>, (void *)args);
        }
    }
    else {
        
        thread_context_t<T> *th_contx = new thread_context_t<T>();
        
       //// th_contx->cur_reads = new std::vector<read_pair_t<T> *> ();
        
        if (numOfpeaks > 0) {
            //debug("init threads init peak_reads");
            th_contx->peak_reads = new std::vector<int>();
            
            for (int j=0; j < numOfpeaks; j++) {
                th_contx->peak_reads->push_back(0);
            }
            th_contx->multi_read_mapTo_pids = new MULTI_READ_MAP ();
            th_contx->peak_uniqReads_startPos = new std::map<int, std::vector<int> > ();
        }
        th_contx->track = new ShortRead();
        
        global_context->thread_contexts.push_back(th_contx);
        
    }
    
    return;
}

//release memory
template <typename T>
void destroy_thread_context(global_context_t<T> * global_context)
{
    
    for(int i=0; i < global_context-> thread_number; i++)
    {
        
        if (global_context->thread_contexts[i]->multi_read_mapTo_pids != nullptr) {
            delete global_context->thread_contexts[i]->multi_read_mapTo_pids;
            delete global_context->thread_contexts[i]->peak_uniqReads_startPos;
        }
        if (global_context->thread_contexts[i]->peak_reads != nullptr) {
            delete global_context->thread_contexts[i]->peak_reads;
        }
        
        
        int size =global_context->thread_contexts[i]->cur_reads.size();
        
        for (int j=0; j< size; j++)
        {
            delete (global_context->thread_contexts[i]->cur_reads[j]);
        }
        global_context->thread_contexts[i]->cur_reads.clear();
        //delete global_context->thread_contexts[i]->cur_reads;
        
        delete global_context->thread_contexts[i]->track;
        
        if(global_context->thread_number >1){
            pthread_spin_destroy(&global_context -> thread_contexts[i]->cur_reads_lock);
        }
        delete global_context -> thread_contexts[i];
    }
    
}

//join threads
template <typename T>
void join_threads(global_context_t<T> * global_context)
{
    int xk1;
    for(xk1=0; xk1<global_context-> thread_number; xk1++)
        pthread_join(global_context -> thread_contexts[xk1]->thread_object, NULL);
}

//merge results
template <typename T>
void merge_results(global_context_t<T> * global_context, ShortRead * track)
{
    long fraglength = 0;
    long taglength = 0;
    int n_frags = 0;
    int n_tags = 0;
    
    for(size_t i=0;i < global_context->thread_number;i++)
    {
        track->append_shortreads(global_context->thread_contexts[i]->track);
        if (global_context->isPE) {
            
            fraglength += global_context->thread_contexts[i]->fraglength;
            n_frags += global_context->thread_contexts[i]->n_frags;
            
        }
        else {
            taglength += global_context->thread_contexts[i]->taglength;
            n_tags += global_context->thread_contexts[i]->n_tags;
        }
    }
    if (global_context->isPE) {
        if (n_frags == 0) {
            error("Error in estimating fragment length : not enough number of fragments.");
            std::exit(1);
        }
        
        fraglength = fraglength/n_frags;
        track->set_fraglength(int(fraglength));
        
        
        //debug("set fraglength " + std::to_string(int(fraglength)));
        
    }
    else {
        if (n_tags == 0) {
            error("Error in estimating fragment length : not enough number of fragments.");
            std::exit(1);
        }
        taglength = taglength/n_tags;
        track->set_taglength(int(taglength));
    }
}

template <typename T>
void merge_results_redistr(global_context_t<T> * global_context, std::vector<double> * peak_reads, MULTI_READ_MAP * multi_read_mapTo_pids, std::map<int, std::vector<int> > * peak_uniqReads_startPos)
{
    
    for(size_t i=0;i < global_context->thread_number;i++)
    {
        //track->append_shortreads(global_context->thread_contexts[i]->track);
        
        std::map<int, std::vector<int>>::iterator itr ;
        for(itr =global_context->thread_contexts[i]->peak_uniqReads_startPos->begin(); itr!=  global_context->thread_contexts[i]->peak_uniqReads_startPos->end(); itr++)
        {
            int pid = itr->first;
            
            if(peak_uniqReads_startPos->find(pid) !=peak_uniqReads_startPos->end() )
            {
                (* peak_uniqReads_startPos)[pid].reserve( (* peak_uniqReads_startPos)[pid].size() + itr->second.size() );                // preallocate memory without erase original data
                (* peak_uniqReads_startPos)[pid].insert( (* peak_uniqReads_startPos)[pid].end(), itr->second.begin(), itr->second.end() );
            }
            else {
                peak_uniqReads_startPos->insert(std::pair<int,std::vector<int> > (pid,itr->second) );
            }
        }
        for (size_t j=0;j < peak_reads->size();j++) {
            
            peak_reads->at(j) += global_context->thread_contexts[i]->peak_reads->at(j);
        }
        
        MULTI_READ_MAP::iterator itr2;
        for (itr2 = global_context->thread_contexts[i]->multi_read_mapTo_pids->begin(); itr2 != global_context->thread_contexts[i]->multi_read_mapTo_pids->end(); itr2 ++)
        {
            multi_read_mapTo_pids->insert(std::pair<int,multi_read_t > (itr2->first,itr2->second));
            
        }
       /* for (itr = global_context->thread_contexts[i]->multi_read_Align_start->begin(); itr != global_context->thread_contexts[i]->multi_read_Align_start->end(); itr ++)
        {
        multi_read_align_starts->insert(std::pair<int,std::vector<int> > (itr->first,itr->second));
        }*/
    }
}

//read alignment files. re-distribute multi-reads among candidate peaks if the last three params are not nullptr
void parse_bam_multi(ShortRead * track,std::string inputFile,int shiftsize,int threadNumber,Candidate_Peaks * peakIdx,std::vector<double> * peak_reads_Prime, MULTI_READ_MAP * multi_read_mapTo_pids,std::map<int, std::vector<int> > * peak_uniqReads_startPos)
{
    
    //global context for multi-threading
    global_context_t<bam1_t> global_context;
    global_context.all_finished = 0;
    global_context.flag = MULTIONLY;
    global_context.shiftsize = shiftsize;
    global_context.peakIdx = peakIdx;
    
    //open SAM/BAM file
    samFile *fp_in = NULL;
    bam1_t * aligned_read = NULL;
    bam_hdr_t *header = NULL;
    
    std::string format = "BAM";
    
    info("start read bam file ...");
    
    fp_in = sam_open(inputFile.c_str(), "rb");
    
    if(NULL == fp_in) {
        
        error("Could not open file " + inputFile);
        std::exit(1);
    }
    else {
        info("open file " + inputFile);
    }
    
    if(hts_get_format(fp_in)->format == sam) {format = "SAM";}
    
    global_context.format = format;
    global_context.isPE = true;
    
    aligned_read = bam_init1();
    
    //Results main_res;
    //read header
    header = sam_hdr_read(fp_in);
    
    if(NULL == header){
        //std::cout << "No header information." << std::endl;
        error("No header information in " + inputFile);
        std::exit(1);
    }
    for(int i=0;i< header->n_targets;i++)
    {
        char * p = header->target_name[i];
        int  l = header->target_len[i];
        std::string chr(p);
        
        global_context.refnames.push_back(chr);
        global_context.reflengths.insert(std::pair<std::string, int>(chr,l));
    }
    
   // debug("read in header information. number of reference sequence is " + std::to_string(global_context.refnames.size()));
    
    std::string prev_read_id = "";
    std::vector<bam1_t *> multi_read1;
    std::vector<bam1_t *> multi_read2;
    
    std::vector<bam1_t *> alignments_per_read;
    
    int lineno = 0;
    int readID = 0;
    
    init_thread<bam1_t>(&global_context, threadNumber);
    
    int current_thread_id = 0;
    
    //read alignment line by line
    while(sam_read1(fp_in, header,aligned_read) > 0)
    {
        //output line number
        lineno +=1 ;
        if (lineno % 1000000 == 0 ) {
            info(std::to_string(lineno) + " alignments.");
        }
        
        //get read name
        std::string qname ="";
        char * name = bam_get_qname(aligned_read);
        if(name !=NULL)
        {
            qname = std::string(name);
        }
        //============================SE reads=============================
        if (! IS_PAIRED(aligned_read))
        {
            global_context.isPE = false;
            //don not count unmapped read
            if (IS_UNMAPPED(aligned_read) ||IS_DUP(aligned_read)|| aligned_read->core.flag & BAM_FQCFAIL )
            {
                bam_destroy1(aligned_read);
                aligned_read =  bam_init1();
                continue;
            }
            //same read as the previous alignment
            if (qname == prev_read_id || prev_read_id == "") {
                alignments_per_read.push_back(aligned_read);
                prev_read_id = qname;
                aligned_read =  bam_init1();
                continue;
            }
            else {
                //new read
                readID += 1;
                read_pair_t<bam1_t> * read_pair = new read_pair_t<bam1_t>();
                
                for(unsigned int i=0;i<alignments_per_read.size();i++){
                    (read_pair->first_reads).push_back(alignments_per_read[i]);
                    
                }
                read_pair->readID = readID;
                
                if(global_context.thread_number >1){
                            
                    thread_context_t<bam1_t> * thread_context = global_context.thread_contexts[current_thread_id];
                            
                    pthread_spin_lock(&thread_context->cur_reads_lock);
                    thread_context->cur_reads.push_back(read_pair);
                    pthread_spin_unlock(&thread_context->cur_reads_lock);
                    current_thread_id++;
                            
                    if(current_thread_id >= global_context.thread_number) current_thread_id = 0;
                }
                else {
                process_aligned_fragment_multi<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                            
                }
                alignments_per_read.clear();
                alignments_per_read.push_back(aligned_read);
                prev_read_id = qname;
                aligned_read = bam_init1();
            }//end new read
        }
        //====================================PE read===============================
        else {//PE read, need to keep tracking of fragment length
            size_t flag_pos = qname.find('/');
            //cur_read_id = aligned_read.qname;
            if (flag_pos != std::string::npos) {qname = qname.substr(0,flag_pos);}
            else {
                if (boost::algorithm::ends_with(qname,".1")|| boost::algorithm::ends_with(qname,".2")) {
                    flag_pos = qname.size() - 2;
                    qname = qname.substr(0,flag_pos);
                }
            }
            //same read as the previous alignment
            if (qname == prev_read_id ){
                if (IS_READ1(aligned_read) ){ multi_read1.push_back(aligned_read);}
                if (IS_READ2(aligned_read) ){multi_read2.push_back(aligned_read);}
                aligned_read = bam_init1();
                continue;
            }
            else {//new read
                //process the previous read
                if (multi_read1.size() >=1 || multi_read2.size() >= 1){
                    readID += 1;
                    read_pair_t<bam1_t> * read_pair = new read_pair_t<bam1_t>();
                    read_pair->readID = readID;
                    
                    for(unsigned int i=0;i< multi_read1.size();i++)
                    {
                        read_pair->first_reads.push_back(multi_read1[i]);
                    }
                    for(unsigned int i=0;i< multi_read2.size();i++)
                    {
                        read_pair->second_reads.push_back(multi_read2[i]);
                    }
                    //multi-threading
                    if(global_context.thread_number > 1){
                        //prepare thread context
                        thread_context_t<bam1_t> * thread_context = global_context.thread_contexts[current_thread_id];
                        
                        pthread_spin_lock(&thread_context->cur_reads_lock);
                        thread_context->cur_reads.push_back(read_pair);
                        pthread_spin_unlock(&thread_context->cur_reads_lock);
                        current_thread_id++;
                        
                        if(current_thread_id >= global_context.thread_number)
                            current_thread_id = 0;
                    }
                    else {
                        process_aligned_fragment_multi<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                        
                    }
                }//end process the previous read
               
                multi_read1.clear();
                multi_read2.clear();
                prev_read_id = qname;
                if (IS_READ1(aligned_read)){
                    multi_read1.push_back(aligned_read);}
                if (IS_READ2(aligned_read)){
                    multi_read2.push_back(aligned_read);}
                
                aligned_read =  bam_init1();
            }//end new read
            
        }//===============end PE reads======================================
    }//end while
    
    //close file
    bam_destroy1(aligned_read);
    bam_hdr_destroy(header);
    sam_close(fp_in);
    
    global_context.all_finished = 1;
    //join threads
    if(global_context.thread_number >1)
    {
        join_threads<bam1_t>(&global_context);
        
    }

    merge_results_redistr<bam1_t>(&global_context,peak_reads_Prime,multi_read_mapTo_pids, peak_uniqReads_startPos);
    
    
    destroy_thread_context<bam1_t>(&global_context);
    
}

//parse alignment files of both uniquely mapped reads and multi-reads.
void parse_bam_both(ShortRead * track,std::string inputFile,int shiftsize,int threadNumber)
{
    
    //global context for multi-threading
    global_context_t<bam1_t> global_context;
    global_context.all_finished = 0;
    global_context.flag = BOTH;
    global_context.peakIdx = nullptr;
    
    
    //open SAM/BAM file
    samFile *fp_in = NULL;
    bam1_t * aligned_read = NULL;
    bam_hdr_t *header = NULL;
    
    std::string format = "BAM";
    
    info("start read bam file ...");
    
    fp_in = sam_open(inputFile.c_str(), "rb");
    
    if(NULL == fp_in) {
        
        error("Could not open file " + inputFile);
        std::exit(1);
    }
    else {
        info("open file " + inputFile);
    }
    
    if(hts_get_format(fp_in)->format == sam) {format = "SAM";}
    
    global_context.format = format;
    global_context.isPE = true;
    
    aligned_read = bam_init1();
    
    //Results main_res;
    //read header
    header = sam_hdr_read(fp_in);
    
    if(NULL == header){
        //std::cout << "No header information." << std::endl;
        error("No header information in " + inputFile);
        std::exit(1);
    }
    for(int i=0;i< header->n_targets;i++)
    {
        char * p = header->target_name[i];
        int  l = header->target_len[i];
        std::string chr(p);
        
        global_context.refnames.push_back(chr);
        global_context.reflengths.insert(std::pair<std::string, int>(chr,l));
    }
    
    
    
    std::string prev_read_id = "";
    std::vector<bam1_t *> multi_read1;
    std::vector<bam1_t *> multi_read2;
    std::vector<bam1_t *> alignments_per_read;
    
    int lineno = 0;
    int readID = 0;
    
    //init threads for parsing alignments
    init_thread<bam1_t>(&global_context, threadNumber);
    int current_thread_id = 0;
    
    //read alignments line by line
    while(sam_read1(fp_in, header,aligned_read) > 0)
    {
        //output line number
        lineno +=1 ;
        if (lineno % 1000000 == 0 ) {
            info(std::to_string(lineno) + " alignments.");
        }
        
        //get read name
        std::string qname ="";
        char * name = bam_get_qname(aligned_read);
        if(name !=NULL)
        {
            qname = std::string(name);
        }
       
        
        //=======================SE reads================
        if (! IS_PAIRED(aligned_read))
        {
            global_context.isPE = false;
            
            //don't count unmapped read
            if (IS_UNMAPPED(aligned_read) ||IS_DUP(aligned_read)|| aligned_read->core.flag & BAM_FQCFAIL )
            {
                bam_destroy1(aligned_read);
                aligned_read =  bam_init1();
                info("a read failed");
                if(IS_UNMAPPED(aligned_read))
                {
                    info("unmapped");
                }
                if(IS_DUP(aligned_read))
                {
                    info("dup");
                }
                if(aligned_read->core.flag & BAM_FQCFAIL)
                {
                    info("FQCfail");
                }
                continue;
            }
            
            if (qname == prev_read_id || prev_read_id == "")
            {//same read as the previous alignment, i.e., multi-reads
                alignments_per_read.push_back(aligned_read);
                prev_read_id = qname;
                aligned_read =  bam_init1();
                continue;
            }
            else {
                //new read
                //process the previous read
                    
                        readID += 1;
                        read_pair_t<bam1_t> * read_pair = new read_pair_t<bam1_t>();
                        
                        for(unsigned int i=0;i<alignments_per_read.size();i++){
                            (read_pair->first_reads).push_back(alignments_per_read[i]);
                        }
                        
                        read_pair->readID = readID;
                        //multi-threading
                        if(global_context.thread_number >1){
                            
                            thread_context_t<bam1_t> * thread_context = global_context.thread_contexts[current_thread_id];
                            
                            pthread_spin_lock(&thread_context->cur_reads_lock);
                            thread_context->cur_reads.push_back(read_pair);
                            pthread_spin_unlock(&thread_context->cur_reads_lock);
                            
                            current_thread_id++;
                            if(current_thread_id >= global_context.thread_number) current_thread_id = 0;
                        }
                        else {
                                process_aligned_fragment<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                        }
                
                //clean the buffer for the new read
                alignments_per_read.clear();
                alignments_per_read.push_back(aligned_read);
                prev_read_id = qname;
                aligned_read = bam_init1();
            }
        }
        //======================PE reads ===============================================
        else {//PE read, need to keep tracking of fragment length
            size_t flag_pos = qname.find('/');
            //cur_read_id = aligned_read.qname;
            if (flag_pos != std::string::npos) {qname = qname.substr(0,flag_pos);}
            else {
                if (boost::algorithm::ends_with(qname,".1")|| boost::algorithm::ends_with(qname,".2")) {
                    flag_pos = qname.size() - 2;
                    qname = qname.substr(0,flag_pos);
                }
            }
            //same as the previous read
            if (qname == prev_read_id ){
                if (IS_READ1(aligned_read) ){ multi_read1.push_back(aligned_read);}
                if (IS_READ2(aligned_read) ){ multi_read2.push_back(aligned_read);}
                aligned_read = bam_init1();
                continue;
            }
            else {//new read
                
                //previous read is a multi-read
                if (multi_read1.size() >=1 || multi_read2.size() >= 1){
                    readID += 1;
                    read_pair_t<bam1_t> * read_pair = new read_pair_t<bam1_t>();
                    read_pair->readID = readID;
                    
                    for(unsigned int i=0;i< multi_read1.size();i++)
                    {
                        read_pair->first_reads.push_back(multi_read1[i]);
                    }
                    for(unsigned int i=0;i< multi_read2.size();i++)
                    {
                        read_pair->second_reads.push_back(multi_read2[i]);
                    }
                    
                    if(global_context.thread_number > 1){
                        
                        thread_context_t<bam1_t> * thread_context = global_context.thread_contexts[current_thread_id];
                        
                        pthread_spin_lock(&thread_context->cur_reads_lock);
                        thread_context->cur_reads.push_back(read_pair);
                        pthread_spin_unlock(&thread_context->cur_reads_lock);
                        current_thread_id++;
                        
                        if(current_thread_id >= global_context.thread_number)
                            current_thread_id = 0;
                    }
                    else {
                    process_aligned_fragment<bam1_t>(&global_context,global_context.thread_contexts[0],read_pair);
                    }
                }//end process previous read
               
                multi_read1.clear();
                multi_read2.clear();
                prev_read_id = qname;
                
                if (IS_READ1(aligned_read)){
                    multi_read1.push_back(aligned_read);}
                if (IS_READ2(aligned_read)){
                    multi_read2.push_back(aligned_read);}
                aligned_read =  bam_init1();
            }//end new read
            
        }//===========================end PE read============================================
    }//end while
    
    //close file
    bam_destroy1(aligned_read);
    bam_hdr_destroy(header);
    sam_close(fp_in);
    
    global_context.all_finished = 1;
    //join threads
    if(global_context.thread_number >1)
    {
        join_threads<bam1_t>(&global_context);
        
    }
 
    track->isPE = global_context.isPE;
    
    merge_results<bam1_t>(&global_context, track);
    //debug("global context reflengths = " + std::to_string(global_context.reflengths.size()));
    track->set_rlengths(global_context.reflengths);
    
    destroy_thread_context<bam1_t>(&global_context);

}

//read in alignments for both uniquely mapped reads and multi-reads
std::pair<ShortRead *, ShortRead *> read_aligmentFile( opt_t &options , int flag)
{
    info("#1 read tags...");
    ShortRead *treat = new ShortRead();
    ShortRead *control = new ShortRead();
    
    if (options.format == "BAM" ) {
        
        parse_bam_both(treat,options.tfile,options.shift,options.threadNum);
        parse_bam_both(control,options.cfile,options.shift,options.threadNum);
        
        if (treat->isPE) {
            options.PE_mode = true;
            options.tsize = treat->fraglength;
        }
        else {
            options.tsize = treat->tsize;
        }
    }
    
    return std::pair<ShortRead *,ShortRead *>(treat, control);
}

//multi-read re-distribution among candidate peaks.
void read_distribution(ShortRead * track, std::string inputFile,Candidate_Peaks * peakIdx,opt_t options,std::vector<double> * peak_reads_Prime, MULTI_READ_MAP * multi_read_mapTo_pids,std::map<int, std::vector<int> > * peak_uniqReads_startPos)
{
    
    if (options.format == "BAM") {
        int shiftsize = options.d/2;
        info("call parse_bam_multi ...");
    parse_bam_multi(track,inputFile,shiftsize,options.threadNum,peakIdx,peak_reads_Prime,multi_read_mapTo_pids,peak_uniqReads_startPos);
    }
    if (options.format == "BED") {
        //int shiftsize = options.d/2;
        info("doesn't support BED file yet.");
        std::exit(0);
       // parse_bed(track,inputFile,shiftsize,options.threadNum, "BED",MULTIONLY,peakIdx,peak_reads_Prime,multi_read_mapTo_pids);
    }
    
    
}


