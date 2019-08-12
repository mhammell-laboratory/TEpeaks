//
//  Fwtrak.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 5/3/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#include <assert.h>
#include <climits>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include "ShortRead.h"
#include "myLog.h"
#include "Pileup.h"
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>


bool pos_comp(pos_t first, pos_t second){
    return first.start < second.start ;
}

ShortRead::ShortRead(){
    //fw = fw;
    //__locations = std::vector<std::map<int,std::vector<int> > >();    // location pairs
    std::map<int,std::vector<int> > plus;
    std::map<int,std::vector<int> > minus;

    __locations.push_back(plus);
    __locations.push_back(minus);
    
    
    std::map<int,std::vector<int> > plus_m;
    std::map<int,std::vector<int> > minus_m;
    
    __locations_m.push_back(plus_m);
    __locations_m.push_back(minus_m);
    
    
    __sorted = false;
    __dup_sorted = false;
    __dup_separated = false;
    max_dup_tag = 1;
    total = 0;           // total tags
    total_uniq = 0;
    total_multi = 0;
    dup_total = 0;           // total tags
    //annotation = anno;   // need to be figured out ??
    //rlengths = {}       // lengths of reference sequences, e.g. each chromosome in a genome
    //buffer_size = buffer_size
    //length = 0
    __destroyed = false;
    fraglength = 0;
    tsize = 0;
    redundant_rate = 0.0;
    isPE = false;

    
}

void ShortRead::set_taglength(int tsize)
{
    this->tsize = tsize;
}

void ShortRead::set_fraglength(int flen)
{
    this->fraglength = flen;
}


void ShortRead::append_shortreads(ShortRead *other)
{
    this->total += other->total;
    this->total_uniq += other->total_uniq;
    this->total_multi += other->total_multi;
    this->length += other->length;
    
    average_template_length = (int)this->length/this->total_uniq;
    
        if (this->isPE) {
            //ADD for (auto p : other->__locations3end[strand]) {
            
            for (auto p : other->__locations_pe) {
                //get the chrom name
                std::string chromosome = other->chromIdToName[p.first];
                
                
                if (this->chromNameToId.find(chromosome) == this->chromNameToId.end()) { //new chromosome
                    
                    int chrID = this->chromNameToId.size();
                    this->chromNameToId.insert(std::pair<std::string,int>(chromosome,chrID));
                    this->chromIdToName.insert(std::pair<int,std::string>(chrID,chromosome));
                    //add new data
                    std::vector<pos_t> poslist;
                    for (auto pos : p.second) {
                        
                        poslist.push_back(pos);
                        
                        //this->__locations_pe[strand][chrID].push_back(pos);
                    }
                    this->__locations_pe.insert(std::pair<int,std::vector<pos_t> >(chrID,poslist));
                    
                }
                else {
                    for (auto pos : p.second) {
                        this->__locations_pe[this->chromNameToId[chromosome]].push_back(pos);
                    }
                }
                
            }
            //ADD for (auto p : other->__locations3end_m[strand]) {
            for (auto p : other->__locations_pe_m) {
                //get the chrom name
                std::string chromosome = other->chromIdToName[p.first];
                
                if (this->chromNameToId.find(chromosome) == this->chromNameToId.end()) { //new chromosome
                    
                    int chrID = this->chromNameToId.size();
                    this->chromNameToId.insert(std::pair<std::string,int>(chromosome,chrID));
                    this->chromIdToName.insert(std::pair<int,std::string>(chrID,chromosome));
                    //add new data
                    std::vector<pos_t> poslist;
                    for (auto pos : p.second) {
                        //this->__locations_pe[strand][chrID].push_back(pos);
                        poslist.push_back(pos);
                        
                    }
                    this->__locations_pe_m.insert(std::pair<int,std::vector<pos_t> >(chrID,poslist));
                }
                else {
                    for (auto pos : p.second) {
                        this->__locations_pe_m[this->chromNameToId[chromosome]].push_back(pos);
                    }
                }
                

            }
        }
        else {
            int strand = 0;
            while (strand < 2) {
                
            for (auto p : other->__locations[strand]) {
                std::string chromosome = other->chromIdToName[p.first];
                
                if (this->chromNameToId.find(chromosome) == this->chromNameToId.end()) { //new chromosome
                    
                    int chrID = this->chromNameToId.size();
                    this->chromNameToId.insert(std::pair<std::string,int>(chromosome,chrID));
                    this->chromIdToName.insert(std::pair<int,std::string>(chrID,chromosome));
                    //add new data
                    std::vector<int> poslist;
                    for (auto pos : p.second) {
                        
                        poslist.push_back(pos);
                        
                        
                        //this->__locations_pe[strand][chrID].push_back(pos);
                    }
                    this->__locations[strand].insert(std::pair<int,std::vector<int> >(chrID,poslist));
                }
                else {
                    for (auto pos : p.second) {
                        this->__locations[strand][this->chromNameToId[chromosome]].push_back(pos);
                    }
                }

            }

            for (auto p : other->__locations_m[strand]) {
                std::string chromosome = other->chromIdToName[p.first];
                
                if (this->chromNameToId.find(chromosome) == this->chromNameToId.end()) { //new chromosome
                    
                    int chrID = this->chromNameToId.size();
                    this->chromNameToId.insert(std::pair<std::string,int>(chromosome,chrID));
                    this->chromIdToName.insert(std::pair<int,std::string>(chrID,chromosome));
                    //add new data
                    std::vector<int> poslist;
                    for (auto pos : p.second) {
                        
                        
                        poslist.push_back(pos);
                        
                        
                        //this->__locations_pe[strand][chrID].push_back(pos);
                    }
                    this->__locations_m[strand].insert(std::pair<int,std::vector<int> >(chrID,poslist));
                }
                else {
                    for (auto pos : p.second) {
                        this->__locations_m[strand][this->chromNameToId[chromosome]].push_back(pos);
                    }
                }
            }
                strand +=1;
        }//end while
            
    }//end if
    
}

void ShortRead::set_rlengths ( std::map<std::string,int> rflengths )
{
       /* """Set reference chromosome lengths dictionary.

        Only the chromosome existing in this fwtrack object will be updated.

        If chromosome in this fwtrack is not covered by given
        rlengths, and it has no associated length, it will be set as
        maximum integer.

        """*/
    
    //set valid_chroms, missed_chroms, extra_chroms;
    //str chrom;

   // valid_chroms = set(__locations.keys()).intersection(rlengths.keys());
    
    for (auto v : chromNameToId)
    {
        
        if (rflengths.find(v.first) == rflengths.end())
        {
            this->rlengths.insert(std::pair<std::string,int> (v.first,INT_MAX));
        }
        else
        {
            
            this->rlengths.insert(std::pair<std::string,int> (v.first,rflengths[v.first]));
        }
       // missed_chroms.push_back(v.first);
    }

    //return true;
}

std::map<std::string,int> ShortRead::get_rlengths ( )
{
        /*"""Get reference chromosome lengths dictionary.

        If self.rlength is empty, create a new dict where the length of
        chromosome will be set as the maximum integer.
        """*/
    if (rlengths.size() ==0){
        for (auto v : chromNameToId)
        {
            rlengths.insert(std::pair<std::string,int>(v.first,INT_MAX));
        }
    }
    return rlengths;
}



std::pair<std::vector<int>,std::vector<double> > ShortRead::pileup_a_chromosome_c ( std::string chrom, std::vector<int> ds, std::vector<double> scale_factor_s, bool uniqOnly,double baseline_value )
{
       /* """pileup a certain chromosome, return [p,v] (end position and value) list.

        This function is for control track. Basically, here is a
        simplified function from FixWidthTrack. We pretend the PE is
        SE data and left read is on plus strand and right read is on
        minus strand.
        
        ds             : tag will be extended to this value to 3' direction,
                         unless directional is False. Can contain multiple extension
                         values. Final pileup will the maximum.
        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        """ */
    
    std::pair<std::vector<int>,std::vector<double> > tmp_pileup, prev_pileup;
    double scale_factor;
    int d, five_shift, three_shift;
    int rlength = rlengths[chrom];
    int chromID = chromNameToId[chrom];

    if (not __sorted){
        sort();
    }

    
    if (ds.size() != scale_factor_s.size()) {
        error("ds and scale_factor_s must have the same length!");
        std::exit(1);
    }

    //prev_pileup = None;
    bool prev_pileup_init_flag=false;
    
    debug("scale_factor_s size  = " + std::to_string(scale_factor_s.size()));
    
    for (size_t i=0; i <  scale_factor_s.size(); i++){
        d = ds[i];
        scale_factor = scale_factor_s[i];
        five_shift = d/2;
        three_shift= d/2;

        std::vector<int> plustags;
        std::vector<int> minustags;
        
        debug(" d = " + std::to_string(d));
        
        for(auto poss : __locations_pe[chromID])
        {
            plustags.push_back(poss.start);
            minustags.push_back(poss.end);
        }
    
        debug("in pileup_a_chromosome_c ");
        tmp_pileup = se_all_in_one_pileup ( plustags, minustags, five_shift, three_shift, rlength, scale_factor, baseline_value );
        
        debug("tmp_pileup size = " + std::to_string(tmp_pileup.first.size()));
        
        
        if (prev_pileup_init_flag){
            
            
            std::pair<std::vector<int>,std::vector<double> >  ret_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup );
            prev_pileup.first.clear();
            debug("prev_pileup size = " + std::to_string(prev_pileup.first.size()));
            prev_pileup.second.clear();
            prev_pileup = ret_pileup;
            debug("prev_pileup size = " + std::to_string(prev_pileup.first.size()));
        }
        else{
            prev_pileup = tmp_pileup;
            prev_pileup_init_flag = true;
        }
    }

    return prev_pileup;
}


std::pair<std::vector<int>,std::vector<double> > ShortRead::pileup_a_chromosome ( std::string chrom, std::vector<int> ds, std::vector<double> scale_factor_s, bool uniqOnly, double baseline_value, bool directional, int end_shift)
{
        /*"""pileup a certain chromosome, return [p,v] (end position and value) list.
        
        ds             : tag will be extended to this value to 3' direction,
                         unless directional is False. Can contain multiple extension
                         values. Final pileup will be the maximum.
        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
        directional    : if False, the strand or direction of tag will be ignored, so that extenstion will be both sides with d/2.
        end_shift      : move cutting ends towards 5->3 direction if value is positive, or towards 3->5 direction if negative. Default is 0 -- no shift at all.

        """ */
    
    //long d;
    long five_shift, three_shift;  // adjustment to 5' end and 3' end positions to make a fragment
    //dict chrlengths = get_rlengths ();
    
    int rlength = rlengths[chrom];
    int chromID = chromNameToId[chrom];
    //object ends;
    std::vector<int> five_shift_s ;
    std::vector<int> three_shift_s ;
    
    if (ds.size() != scale_factor_s.size()) {
        error("ds and scale_factor_s must have the same length!");
        std::exit(1);
    }
    
    //adjust extension length according to 'directional' and 'halfextension' setting.
    for (auto d : ds){
        debug("d  = " + std::to_string(d));
        if (directional){
                // only extend to 3' side
            five_shift_s.push_back(  - end_shift );
            three_shift_s.push_back( end_shift + d);
        }
        else{
                // both sides
            five_shift_s.push_back( d/2 - end_shift );
            three_shift_s.push_back( end_shift + d - d/2);
        }
    }
    
    std::pair<std::vector<int>,std::vector<double> >  prev_pileup;
    
    bool prev_pileup_init_flag = false;
    //for( size_t i =0; i < ds.size() ; i++)

    for( size_t i =0; i <  ds.size(); i++)
    {
        five_shift = five_shift_s[i];
        three_shift = three_shift_s[i];
        double scale_factor = scale_factor_s[i];
        //SE
        std::pair<std::vector<int>,std::vector<double> > tmp_pileup;
        
        
        if (uniqOnly) {
            
            tmp_pileup = se_all_in_one_pileup ( __locations[0][chromID], __locations[1][chromID], five_shift, three_shift, rlength, scale_factor, baseline_value );
            debug("tmp_pileup size = " + std::to_string(tmp_pileup.first.size()));
            
        }
        else {//include multi-reads
            std::vector<int> merge_list1; //merge unique reads and multi reads, strand=0
            std::vector<int> merge_list2; //merge unique reads and multi reads, strand = 1

            std::merge(__locations[0][chromID].begin(),__locations[0][chromID].end(),__locations_m[0][chromID].begin(),__locations_m[0][chromID].end(),std::back_inserter(merge_list1));
            
            std::merge(__locations[1][chromID].begin(),__locations[1][chromID].end(),__locations_m[1][chromID].begin(),__locations_m[1][chromID].end(),std::back_inserter(merge_list2));
            
            tmp_pileup = se_all_in_one_pileup ( merge_list1, merge_list2, five_shift, three_shift, rlength, scale_factor, baseline_value );
        }

        if (prev_pileup_init_flag ){
           std::pair<std::vector<int>,std::vector<double> >  ret_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup );
            prev_pileup.first.clear();
            prev_pileup.second.clear();
            prev_pileup = ret_pileup;
            
            debug("prev_pileup size = " + std::to_string(prev_pileup.first.size()));
        }
        else{
            prev_pileup = tmp_pileup;
            prev_pileup_init_flag = true;
        }
    }
    return prev_pileup;
}

std::vector<std::string> ShortRead::get_chrom_names(){
    std::vector<std::string> ret;
    
    for (auto p : this->chromNameToId) {
        ret.push_back(p.first);
    }

   return ret;
}

std::pair<std::vector<int>,std::vector<double> > ShortRead::pileup_a_chromosome_pe ( std::string chrom, std::vector<double> scale_factor_s, bool uniqOnly, double baseline_value )
{
    /*"""pileup a certain chromosome, return [p,v] (end position and value) list.
        
        scale_factor_s  : linearly scale the pileup value applied to each d in ds. The list should have the same length as ds.
        baseline_value : a value to be filled for missing values, and will be the minimum pileup.
    */
    
    std::pair<std::vector<int>,std::vector<double> > tmp_pileup, prev_pileup;
    double scale_factor;

    int chromID = chromNameToId[chrom];
    bool prev_pileup_init_flag = false;
    
    
    for (size_t i =0; i < scale_factor_s.size(); i++)
    {
        scale_factor = scale_factor_s[i];

        tmp_pileup = quick_pileup ( __locations_pe[chromID], scale_factor, baseline_value );
        
        
        if (prev_pileup_init_flag){
            prev_pileup = max_over_two_pv_array ( prev_pileup, tmp_pileup );
            
        }
        else{
            prev_pileup = tmp_pileup;
            prev_pileup_init_flag = true;
        }
    }
    return prev_pileup;
}


int ShortRead::get_num_of_chroms()
{
    return chromNameToId.size();
}

//return number of unique reads in chromosome chr, used by PeakModel, for SE only
int ShortRead::numOfUniqTags_by_chr(int chrID, int strand)
{
    if (chromIdToName.find(chrID) != chromIdToName.end()) {
        return __locations[strand][chrID].size();
    }
    else {
        error( "Chromosome " + chromIdToName[chrID]  + " does not exist in the alignment file.");
        std::exit(1);
    }
}

//return locations of unique reads, used by PeakModel, for SE only
std::vector<int> ShortRead::get_locations_by_chr_uniq(int chrID, int strand)
{
    if (chromIdToName.find(chrID) != chromIdToName.end()) {
        return __locations[strand][chrID];
    }
    else {
        error("chromosome " + chromIdToName[chrID] + " does not exist.");
        std::exit(1);
    }
}

ShortRead::~ShortRead(){
    
}
//Single-end reads: add weight to alignment, one position could have more than one read
void ShortRead::add_loc ( std::string chromosome, int fiveendpos, int strand , double w)
{
        /*"""Add a location to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        strand     -- 0: plus, 1: minus
        """*/

    int chrID = -1;
    
    if (chromNameToId.find(chromosome) == chromNameToId.end()) { //new chromosome
        
        chrID = chromNameToId.size();
        chromNameToId.insert(std::pair<std::string,int>(chromosome,chrID));
        chromIdToName.insert(std::pair<int,std::string>(chrID,chromosome));
        
    }
    else {
        chrID = chromNameToId[chromosome];
    }
    
    if (w == 1) { //unique reads

        if (__locations[strand].find(chrID) == __locations[strand].end())
        {
            std::vector<int> pos;
            pos.push_back(fiveendpos);
            
            __locations[strand].insert(std::pair<int,std::vector<int> > (chrID,pos));

        }
        else
        {
            __locations[strand][chrID].push_back(fiveendpos);
        }
        
    }
    else { //multi reads
        if (__locations_m[strand].find(chrID) == __locations_m[strand].end())
        {
            std::vector<int> pos;
            pos.push_back(fiveendpos);
            
            __locations_m[strand].insert(std::pair<int,std::vector<int> >(chrID, pos));
            
        }
        else
        {
            __locations_m[strand][chrID].push_back(fiveendpos);
        }
        
        
    }
}


//Paired-end reads : add weight to alignment, one position could have more than one read,
void ShortRead::add_loc_pe ( std::string chromosome, int fiveendpos, int threeendpos, double w )
{
    //5 end pos, 3 end pos, and weight
    int chrID = -1;
    
    if (chromNameToId.find(chromosome) == chromNameToId.end()) { //new chromosome
        
        
        chrID = chromNameToId.size();
        chromNameToId.insert(std::pair<std::string,int>(chromosome,chrID));
        
        chromIdToName.insert(std::pair<int,std::string>(chrID,chromosome));
        
    }
    else {
        chrID = chromNameToId[chromosome];
    }
    
   
    pos_t pos;
    if(threeendpos < fiveendpos)
    {
        pos.start = threeendpos;
        pos.end = fiveendpos;
    }
    else {
        pos.start = fiveendpos;
        pos.end = threeendpos;
    }
    
    if (w == 1.0) { //unique reads
        
        if (__locations_pe.find(chrID) == __locations_pe.end()){
    
            std::vector<pos_t> poslist;
            poslist.push_back(pos);
        
            __locations_pe.insert(std::pair<int,std::vector<pos_t> >(chrID,poslist));
            
        }
        else
        {
            __locations_pe[chrID].push_back(pos);
            
            
        }
        this->length += pos.end - pos.start;
        
        if(this->total_uniq != 0 ){
            this->average_template_length = (int) this->length/this->total_uniq;
        }
    }
    else {
        //multi-reads
        if (__locations_pe_m.find(chrID) == __locations_pe_m.end()){
            
            std::vector<pos_t> poslist;
            
            poslist.push_back(pos);
            
            __locations_pe_m.insert(std::pair<int,std::vector<pos_t> >(chrID,poslist));
            
        }
        else
        {
            
            __locations_pe_m[chrID].push_back(pos);
        }
        
    }
}

void ShortRead::clearMultiReads()
{
    __locations_m.clear();
    __locations_pe_m.clear();
}

void ShortRead::inc_uniq_read()
{
    this->total_uniq += 1;
    this->total += 1;
}
void ShortRead::inc_multi_read()
{
    this->total_multi += 1;
    this->total +=1;
}

int ShortRead::get_total()
{
    return this->total_multi + this->total_uniq + 1;
}
int ShortRead::get_total_uniq()
{
    return this->total_uniq + 1;
}

void ShortRead::clean_m()
{
    if (isPE) {
        __locations_pe_m.clear();
    }
    else {
        __locations_m.clear();
    }
   
}

void ShortRead::set_redundant_rate(double val)
{
    this->redundant_rate = val;
    this->total_multi = (int) this->total_multi * (1.0 - val);
}

void ShortRead::sort()
{
    
    if (isPE) {
        for (auto p : __locations_pe) {
            std::sort(p.second.begin(),p.second.end(),pos_comp);
        }
        for (auto p : __locations_pe_m) {
            std::sort(p.second.begin(),p.second.end(),pos_comp);
        }
        return;
    }
    else {
        for(auto chr_nameID : chromNameToId)
        {
            std::sort(__locations[0][chr_nameID.second].begin(),__locations[0][chr_nameID.second].end());
            std::sort(__locations[1][chr_nameID.second].begin(),__locations[1][chr_nameID.second].end());
            
            std::sort(__locations_m[0][chr_nameID.second].begin(),__locations_m[0][chr_nameID.second].end());
            std::sort(__locations_m[1][chr_nameID.second].begin(),__locations_m[1][chr_nameID.second].end());
        }
        
        return ;
    }

}

void ShortRead::separate_dups( int maxint )
{
     // """Separate the duplicated reads into a different track stored at dup

    this->max_dup_tag = maxint;
    
    if (not __sorted){
        this->sort();
    }

    //__dup_pointer = copy(self.__pointer);
    this->dup_total = 0;
    this->total_uniq = 0;
    this->length = 0;
    
    for(auto chr_nameID : chromNameToId)
    {
        //for each chromosome.
        std::string chrom = chr_nameID.first;
        int i_new = 0;
        int i_dup = 0;
        size_t i_old=1;
        
        // plus strand
        int size = __locations[0][chr_nameID.second].size(); //number of distinctive position
        
        if (size > 1){
           
            // first item
            i_new += 1;
            int current_loc = __locations[0][chr_nameID.second][0];
            int n = 1;
            
            while(i_old < __locations[0][chr_nameID.second].size())
            {
                int p = __locations[0][chr_nameID.second][ i_old ];
               
                if (p == current_loc)
                {
                    n += 1;
                }
                else{
                    current_loc = p;
                    n = 1;
                }
                if (n > maxint)
                {
                    //dup_plus.push_back(p);
                    i_dup += 1;
                    __locations[0][chr_nameID.second].erase(__locations[0][chr_nameID.second].begin()+i_old);
                }
                else {
                     i_old += 1;
                }
                
            }

        }
        this->total_uniq += __locations[0][chr_nameID.second].size();
        this->dup_total += i_dup;

        //- strand
        i_dup = 0;
        
        i_old = 1;
        if (__locations[1][chr_nameID.second].size() > 1){
            
           // new_minus.push_back(minus[ i_new ]);// # first item
           // i_new += 1;
            int current_loc = __locations[1][chr_nameID.second][0];
            int n = 1;
            
            while ( i_old < __locations[1][chr_nameID.second].size())
            {
                int p = __locations[1][chr_nameID.second][ i_old ];
                
                if (p == current_loc){
                    n += 1;
                }
                else{
                    current_loc = p;
                    n = 1;
                }
                if (n > maxint){
                    
                    __locations[1][chr_nameID.second].erase(__locations[1][chr_nameID.second].begin()+i_old);
                    i_dup += 1;
                }
                else {
                    i_old += 1;
                }
                
            }

        }
        
        this->total_uniq +=  __locations[1][chr_nameID.second].size();
        this->dup_total +=  i_dup;
    }
    __dup_separated = true;
    
    
}
/*
void ShortRead::separate_dups_pe (  int maxint )
{
    //"""Filter the duplicated reads. Run it right after you add all data into this object.

    int n, start, end, size;
    
    this->max_dup_tag = maxint;
    
    if (not __sorted){
        this->sort();
    }
    
    //__dup_pointer = copy(self.__pointer)
    this->dup_total = 0;
    this->total_uniq = 0;
    this->length = 0;
    this->average_template_length = 0.0;

    for( auto chr_nameID : chromNameToId)
    { //# for each chromosome
        int i_new = 0;
        int i_dup = 0;
        int i_old = 1 ;
        
        std::vector<pot_t> plus = __locations_pe[0][chr_nameID.second];
        int size = plus.size();
        
        if (size > 1){
            n = 1;
            pos_t  cur_pos = plus[0];

            this->length += cur_pos.end - cur_pos.start;
            
            while (i_old < plus.size())
            {
                pos_t pos = plus[i_old];
                
                bool all_same = false ;
                if((pos.start == cur_pos.start) and (pos.end == cur_pos.end))
                {
                    all_same = true;
                }
                if (all_same){
                    n += 1;
                }
                else{
                    cur_pos = pos;
                    n = 1;
                }
                if (n > maxint){

                    i_dup += 1;
                    plus.erase(plus.begin()+i_old);
                }
                else{
                    
                    this->length += pos.end - pos.start;
                    i_old += 1;
                }
            }
            
            this->dup_total += i_dup;
        }
        this->total_uniq += plus.size();
        
        //- strand
        i_dup = 0;
        
        i_old = 1;
        std::vector<pot_t> minus = __locations_pe[1][chr_nameID.second];
        int size = minus.size();
        if (size > 1){
            
            // new_minus.push_back(minus[ i_new ]);// # first item
            // i_new += 1;
            int current_loc = __locations[1][chr_nameID.second][0];
            int n = 1;
            
            while ( i_old < __locations[1][chr_nameID.second].size())
            {
                
                pos_t pos = minus[i_old];
                bool all_same = false ;
                if((pos.start == cur_pos.start) and (pos.end == cur_pos.end))
                {
                    all_same = true;
                }
                if (all_same){
                    n += 1;
                }
                else{
                    cur_pos = pos;
                    n = 1;
                }
                
                if (n > maxint){
                    
                    minus.erase(minus.begin()+i_old);
                    
                    i_dup += 1;
                }
                else {
                    i_old += 1;
                }
                
            }
            
        }
        
        this->total_uniq +=  __locations[1][chr_nameID.second].size();
        this->dup_total +=  i_dup;
        
    }
            //self.__dup_locations[k] = dup_locs
    this->average_template_length = 1.0 * this->length  / this->total_uniq;
    
    
} */

void ShortRead::sample_num (int samplesize, int seed )
{
    //"""Sample the tags for a given percentage. Warning: the current object is changed!

    double percent = double(samplesize)/get_total();
    this->sample_percent ( percent, seed );
}

void ShortRead::sample_percent (double percent, int seed )
{
    //"""Sample the tags for a given percentage. Warning: the current object is changed!
    
    this->total_uniq = 0;
    this->length = 0;
    
    if (seed >= 0){
        std::srand(seed);
    }
    else {
        std::srand ( unsigned ( 0 ) );
    }

    for (auto chr_nameID : chromNameToId){
        //for each chromosome.
        //This loop body is too big, I may need to split code later...

        if(this->isPE){
            
            int num = (int)(__locations_pe[chr_nameID.second].size() * percent);
            std::random_shuffle(__locations_pe[chr_nameID.second].begin(),__locations_pe[chr_nameID.second].end());
            __locations_pe[chr_nameID.second].resize(num);
            std::sort(__locations_pe[chr_nameID.second].begin(),__locations_pe[chr_nameID.second].end(),pos_comp);
            
            this->total_uniq += num;
            for (auto p : __locations_pe[chr_nameID.second]){
                this->length += p.end - p.start;
            }
            
        }
        else {
            int num1 = (int)(__locations[0][chr_nameID.second].size()* percent);
            std::random_shuffle(__locations[0][chr_nameID.second].begin(),__locations[0][chr_nameID.second].end());
            __locations[0][chr_nameID.second].resize(num1);
            std::sort(__locations[0][chr_nameID.second].begin(),__locations[0][chr_nameID.second].end());
            
            int num2 = (int)(__locations[1][chr_nameID.second].size()* percent);
            std::random_shuffle(__locations[1][chr_nameID.second].begin(),__locations[1][chr_nameID.second].end());
            __locations[1][chr_nameID.second].resize(num2);
            std::sort(__locations[1][chr_nameID.second].begin(),__locations[1][chr_nameID.second].end());
            
            this->total_uniq = num1 + num2;
        }
    }
    
    this->average_template_length = (int)this->length/this->total_uniq;

}



