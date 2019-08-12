//
//  Candidate_Peaks.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/22/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#include "Candidate_Peaks.h"

//#include "IntervalTree.h"
#include <math.h>
#include <fstream>
#include <sstream>
//#include <regex>
#include "stdlib.h"
#include <algorithm>
#include <iostream>


Candidate_Peaks::Candidate_Peaks(std::string peakfilename){
    
    read_peaks(peakfilename);
    
}

int Candidate_Peaks::get_numofpeaks()
{
    return peak_reads.size();
}
//Reading & processing peak files
void Candidate_Peaks::read_peaks(std::string peak_filename)
{
    //chrom -> (peak_id -> start"\t"end)
    std::map<std::string, std::vector<std::pair<int, std::string> > > temp_map ;
    
    std::map<std::string, std::vector<std::pair<int, std::string> > >::iterator tmp_itr;

    int peak_id = 0 ;
    int line_no = 0;
    //int summit = -1;
    int start = -1;
    int end = -1;
   // int length = -1;
    double log10pval = 0.0;
    //double fe = 0;
    //double fdr = -1;
    double tags = -1;
    
    std::ifstream input;
    
    try{
        input.open (peak_filename, std::ifstream::in);

        while(! input.eof()){
            
            std::string line,chrom,start_ss,end_ss,length_ss,summit_ss,tags_ss,log10pval_ss,fe_ss,fdr_ss;
            
            std::stringstream ss;
            
            if (! std::getline(input,line)){
                break;
            }
            line_no ++;
            
            if (line == "\n" || line == "" || !line.compare(0,1,"#") || line_no == 1) {
                continue;
            }
            ss << line;
            
            //std::cout << line << std::endl;
            
            try{
            std::getline(ss,chrom,'\t');
            
            if (chrom == "chrom") {
                continue;
            }
            std::string peakid;
            std::getline(ss,start_ss,'\t');
            std::getline(ss,end_ss,'\t');
            std::getline(ss,peakid,'\t');
            //std::getline(ss,length_ss,'\t');
            //std::getline(ss,summit_ss,'\t');
            std::getline(ss,tags_ss,'\t');
            std::getline(ss,log10pval_ss,'\t');
            //std::getline(ss,fe_ss,'\t');
            //std::getline(ss,fdr_ss,'\t');

                start = std::stoi(start_ss);
                end = std::stoi(end_ss);
                //summit = std::stoi(summit_ss);
                tags = std::stod(tags_ss);
                log10pval = std::stod(log10pval_ss);
                //fe = std::stod(fe_ss);
                //fdr = std::stod(fdr_ss);
                //length = std::stoi(length_ss);
                
                
                peak_reads.push_back(tags);
                peak_start.push_back(start);
                peak_end.push_back(end);
                peak_length.push_back(end - start + 1);
                peak_chrom.push_back(chrom);
                
                //peak_input_reads.push_back(0);
                peak_pval.push_back(log10pval);
                //peak_fe.push_back(fe);
                //peak_summit.push_back(summit);
                
                
                tmp_itr = temp_map.find(chrom);
                
                if (tmp_itr != temp_map.end()) {
                    temp_map[chrom].push_back(std::pair<int,std::string>(peak_id,start_ss+"\t"+end_ss));
                }
                else{
                    std::vector<std::pair<int,std::string> > tmp_plist;
                    tmp_plist.push_back(std::pair<int,std::string> (peak_id,start_ss+"\t"+end_ss));
                    temp_map.insert(std::pair<std::string,std::vector<std::pair<int,std::string> > > (chrom,tmp_plist));
                    
                }
                peak_id += 1;
                
            }
            catch (const std::invalid_argument& ia) {
                std::cerr << "Invalid argument: " << ia.what() << '\n';
                std::exit(1);
            }

            
        }
        
        input.close();
        
    }
    catch(std::ifstream::failure e){
        std::cout << "error in reading file " << peak_filename << std::endl;
    }
    
    if (peak_id == 0 ){
        std::cout << "Warning: No peaks found in peak file." << std::endl;
        std::exit(1);
    }
    
    build_tree(temp_map);

    
}


double Candidate_Peaks::get_count(int g)
{
    if ((size_t)g < peak_reads.size() && g > -1) {
        return peak_reads[g];
    }
    else {
        return -1.0;
    }
}

int Candidate_Peaks::get_length(int g)
{
    if ((size_t)g < peak_length.size() && g > -1) {
        return peak_length[g];
    }
    else {
        return -1;
    }
}

std::string Candidate_Peaks::get_chrom(int g)
{
    if ((size_t)g < peak_chrom.size() && g > -1) {
        return peak_chrom[g];
    }
    else {
        return "";
    }
}

int Candidate_Peaks::get_start(int g)
{
    if ((size_t)g < peak_start.size() && g > -1) {
        return peak_start[g];
    }
    else {
        return -1;
    }
}

int Candidate_Peaks::get_end(int g)
{
    if ((size_t)g < peak_end.size() && g > -1) {
        return peak_end[g];
    }
    else {
        return -1;
    }
}


double Candidate_Peaks::get_pval(int g)
{
    if ((size_t) g < peak_pval.size() && g > -1) {
        return peak_pval[g];
    }
    else {
        return -1;
    }
}
/*int Candidate_Peaks::get_summit(int g)
{
    if (g < peak_summit.size()) {
        return peak_summit[g];
    }
    else {
        return -1;
    }
}
*/
double Candidate_Peaks::get_fe(int g)
{
    if ((size_t) g < peak_fe.size() && g > -1) {
        return peak_fe[g];
    }
    else {
        return -1;
    }
}


/*bool Candidate_Peaks::decrease_cnt(int g)
{
    if (g < peak_reads.size()) {
        peak_reads[g] -= 1;
        return true;
    }
    else {
        return false;
    }
}
*/


Candidate_Peaks::~Candidate_Peaks(){
    chrom_itvTree_Dict_itr it;
    
    for (it=idx_peaks.begin(); it != idx_peaks.end(); it++) {
        std::map<int,IntervalTree *> tmp = it->second;
        std::map<int,IntervalTree *>::iterator tmp_itr ;
        //std::cout << tmp.size() << std::endl;
        for (tmp_itr = tmp.begin(); tmp_itr != tmp.end(); tmp_itr ++) {
            //std::cout << "delete tree " << std::endl;
            //std::cout << tmp_itr->first << std::endl;
            //std::cout << tmp_itr->second->center << std::endl;
            delete tmp_itr->second;
        }
        
    }
    
    
}

void Candidate_Peaks::build_tree(std::map<std::string, std::vector<std::pair<int,std::string> > > peak_chrom_maps)
{
    std::vector<Interval> itemlist;
    std::vector<Interval> sublist;
    std::map<std::string, std::vector<std::pair<int,std::string> > >::iterator it;
    std::vector<std::pair<int,std::string> > tmp;
    
    //gene_exon_Dict_It exon_str_itr;
    //Gene tmp_gene;
    
    std::string chr,pid;//,start_ss,end_ss;
    std::stringstream ss;
    
    //long st, end;
    int  cur_bin_id,start_bin_id, js, je, k; //, buket_size;
    
    for (it = peak_chrom_maps.begin(); it != peak_chrom_maps.end(); it++) {
        chr = it->first;
        tmp = it->second;
        itemlist.clear();
            
        for (size_t i=0; i< tmp.size(); i++) {
                std::string start_ss,end_ss;
                std::stringstream ss,id;
                int start, end;
                id << tmp[i].first;
                ss << tmp[i].second;
                id >> pid;
                
                std::getline(ss,start_ss,'\t');
                std::getline(ss,end_ss,'\t');
                
                start = std::stoi(start_ss);
                end = std::stoi(end_ss);
                
                itemlist.push_back(Interval(std::stoi(pid),-1,start,end,0));
                
            }
        
        std::sort(itemlist.begin(),itemlist.end(),itv_comp) ; //key=operator.attrgetter('start'));
        //quick_sort(itemlist,0,itemlist.size()-1);
        start_bin_id = itemlist[0].start/BIN_SIZE;
        js = 0 ;
        je = 0 ;
        k = 0 ;
    
        for (size_t i=0; i < itemlist.size(); i++) {
            cur_bin_id = itemlist[i].start/BIN_SIZE;
            //std::cout << cur_bin_id   << std::endl;
            
            if (cur_bin_id == start_bin_id) {
                je += 1;
            }
            else {
                //buket_size = (int)sqrt(je - js) + 1;
                //std::cout << buket_size << std::endl;
                
                sublist = std::vector<Interval>(itemlist.begin()+js,itemlist.begin() + je);

                idx_peaks[chr][start_bin_id] = new IntervalTree(sublist);
                // std::cout << start_bin_id << " built one tree."  << std::endl;
                k ++;
                start_bin_id = cur_bin_id ;
                js = je ;
                je ++;
            }
        }
        
        if (js != je) {
            //buket_size = (int) sqrt(je - js) + 1;
            sublist = std::vector<Interval>(itemlist.begin()+js,itemlist.begin() + je);

            idx_peaks[chr][start_bin_id] = new IntervalTree(sublist);
            //print("tree depth = " + str(cds_exon_idx_plus[chr][start_bin_id].get_depth())+ "\n")
            k+=1;
        }
        
    }
    
}

int Candidate_Peaks::get_ovp_peaks(std::string chrom,std::vector<std::pair<int,int> > itv_list, int shiftsize,std::string strand)
{
    //std::vector<std::string> genes;
    std::vector<int> fs ;
   // std::vector<Interval> tmp ;
    chrom_itvTree_Dict::iterator chrom_it;
    std::map<int,IntervalTree *>::iterator bin_iter;
    
    std::string type="" ;
    int bin_id_s, bin_id_e;
    
    for (size_t i=0; i < itv_list.size(); i ++) {
        int start = itv_list[i].first;
        int end = itv_list[i].second;
        
        if (strand == "+") {
            start += shiftsize;
            end += shiftsize;
        }
        if (strand == "-") {
            start = 0 ? start < shiftsize : start - shiftsize;
            end = 0 ? end < shiftsize : end - shiftsize;
        }
        bin_id_s = start/BIN_SIZE;
        bin_id_e = end/BIN_SIZE;
        
        
        chrom_it = idx_peaks.find(chrom);
        if (chrom_it !=  idx_peaks.end()) {
                
                bin_iter = idx_peaks[chrom].find(bin_id_s);
                if (bin_iter != idx_peaks[chrom].end()  )
                {
                    //std::cout << "start to search tree" << std::endl;
                    fs = (*idx_peaks[chrom][bin_id_s]).find_gene(itv_list[i].first,itv_list[i].second);
                    
                    if (bin_id_s != bin_id_e) {
                        bin_iter = idx_peaks[chrom].find(bin_id_e);
                        
                        if (bin_iter != idx_peaks[chrom].end()  ){
                        std::vector<int>  tmp = (*idx_peaks[chrom][bin_id_e]).find_gene(itv_list[i].first,itv_list[i].second);
                            
                        fs.insert(fs.end(),tmp.begin(),tmp.end());
                        }
                    }
                    //std::cout << fs.size() << std::endl;
                }
            }
        
    }
    int pid = -1;
    
    for (size_t i=0; i < fs.size(); i++) {
        pid = fs[i];
        return pid;
    }
    return pid;

}

/*
int main() {
    std::string filename = "test_peak.bed";
    
    Candidate_Peaks peakIdx(filename);
    int numOfpeaks = peakIdx.get_numofpeaks();
     
 //std::vector<chr_ITV> itv_list ;
 //chr_ITV exp ;
    std::string chrom = "chr1";
 int start = 11151;
 int end = 11251;
 //itv_list.push_back(exp);
    std::vector<std::pair<int,int> > itv_list;
    itv_list.push_back(std::pair<int,int> (start,end));

    std::cout << numOfpeaks << std::endl;
    for (int i=0; i< numOfpeaks; i++) {
        std::cout << peakIdx.get_chrom(i) << "\t" << peakIdx.get_start(i) << "\t" << peakIdx.get_end(i) << "\t" << peakIdx.get_count(i) << "\t" << peakIdx.get_fe(i) << "\t" << peakIdx.get_pval(i) << std::endl;
    }
    
    int pid = peakIdx.get_ovp_peaks(chrom,itv_list, 50,"+");
    
    std::cout << pid << "\t" << peakIdx.get_start(pid) << std::endl;

 
 }*/
