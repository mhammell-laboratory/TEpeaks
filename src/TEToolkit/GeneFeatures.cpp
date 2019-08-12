//
//  GeneFeatures.cpp
//  BAMQC-0.5
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include "GeneFeatures.h"

#include <math.h>
#include <fstream>
#include <sstream>
//#include <regex>
#include "stdlib.h"
#include <algorithm>
#include <iostream>
#include <math.h>
#include <iterator>


//#include <boost/tokenizer.hpp>

//template <class T1, class T2, class Pred = std::less<T2> >
//struct sort_pair_second {
//(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
//        Pred p;
//        return p(left.second, right.second);
//    }
//};
//bool sort_pair_second(std::pair<int,int> first, std::pair<int,int> second)
//{
//   return first.second > second.second ;
//}



int pivot(std::vector<Interval> &intervals, int first, int last)
{
    int  p = first;
    int pivotElement = intervals[first].start;
    
    
    for(int i = first+1 ; i <= last ; i++)
    {
        /* If you want to sort the list in the other order, change "<=" to ">" */
        if(intervals[i].start <= pivotElement)
        {
            std::swap(intervals[i],intervals[p]);
            p++;
            
        }
    }
    
    return p;
}

void quick_sort(std::vector<Interval> &intervals, int first, int last){
    
    int pivotElement;
    
    if(first < last)
    {
        pivotElement = pivot(intervals, first, last);
        quick_sort(intervals, first, pivotElement-1);
        quick_sort(intervals, pivotElement+1, last);
    }
    
}


bool reverse_ord_func (int i,int j) { return (j<i); }


Gene::Gene(std::string gid,  std::string ss){
    id = gid;
    strand = ss ;
    min_start =  std::numeric_limits<int>::max();
    max_stop =  std::numeric_limits<int>::min();
    gene_actual_len = 0;
    stop_codon_st = -1;
    stop_codon_end = -1;
}

Gene::~Gene(){}

void Gene::add_cds(int st, int end){
    std::pair<int,int> cds_interval (st,end);
    cds.push_back(cds_interval);
}

void Gene::add_exons(int st, int end){
    if (this->min_start > st) {
        this->min_start = st;
    }
    if (this->max_stop<end) {
        this->max_stop = end;
    }
    std::pair<int,int> exon_interval (st,end);
    gene_actual_len += (end - st+1);
    exons.push_back(exon_interval);
}
void Gene::set_stop_codon(int st,int end)
{
    stop_codon_st = st;
    stop_codon_end = end;
}


void Gene::get_others(){
    std::vector<std::pair<int,int> > left_cds;
    std::vector<std::pair<int,int> > left_exons;
    //int idx[exons.size()];
    size_t i; //,j;
    int itgUp1k_st, itgUp1k_end,itgDn1k_st,itgDn1k_end ;

    sort(exons.begin(),exons.end());
    sort(cds.begin(),cds.end());
    

    for (i = 1; i < exons.size(); i++) {
        intron.push_back(std::make_pair(exons[i-1].second +1, exons[i].first -1));
    }

    if(strand == "+") {
        itgUp1k_st = std::max(int(0),exons[0].first-1000);
        itgUp1k_end = std::max(int(0),exons[0].first-1 );
        itgDn1k_st = exons[exons.size()-1].second + 1;
        itgDn1k_end = exons[exons.size()-1].second + 1000 ;
        itg1k.push_back(std::make_pair(itgUp1k_st,itgUp1k_end)) ;
        itg1k.push_back(std::make_pair(itgDn1k_st,itgDn1k_end)) ;
    
        if (stop_codon_st == -1) {
            utr5 = exons;
        }
        else {
        cds[cds.size()-1].second = stop_codon_end;  
        for( i=0;i < exons.size(); i++) {
                
                if (exons[i].second < cds[0].first) { utr5.push_back(exons[i]); }
            
                if (exons[i].first < cds[0].first && exons[i].second > cds[0].first ) { utr5.push_back(std::make_pair(exons[i].first,cds[0].first - 1)); }
                
                if (exons[i].first <= stop_codon_st && exons[i].second > stop_codon_end )
                { utr3.push_back(std::make_pair(stop_codon_end + 1,exons[i].second)); }
                if (exons[i].first > stop_codon_end ) { utr3.push_back(exons[i]) ; }
            
        }
        }

    }
    else {
        itgDn1k_st = std::max(int(0),exons[0].first - 1000);
        itgDn1k_end = std::max(int(0),exons[0].first -1);
        itgUp1k_st = exons[exons.size()-1].second + 1 ;
        itgUp1k_end = exons[exons.size()-1].second + 1000 ;
        itg1k.push_back(std::make_pair(itgUp1k_st,itgUp1k_end));
        itg1k.push_back(std::make_pair(itgDn1k_st,itgDn1k_end));

        if (stop_codon_st == -1 ) { utr3 = exons; }
        else {
        //if (left_cds.size() == 1) {
           cds[0].first = stop_codon_st;  
           for( i=0;i < exons.size(); i++) {
                
                if (exons[i].second < stop_codon_st) { utr3.push_back(exons[i]); }
            
                if (exons[i].first < stop_codon_st && exons[i].second >= stop_codon_end )
                { utr3.push_back(std::make_pair(exons[i].first,stop_codon_st - 1)); }
                
                if (exons[i].first <= cds[cds.size()-1].first && exons[i].second > cds[cds.size()-1].second )
                { utr5.push_back(std::make_pair(cds[cds.size()-1].second + 1,exons[i].second)); }
            
                if (exons[i].first > cds[cds.size()-1].second ) { utr5.push_back(exons[i]) ; }
           }
        }
        
    }

}



GeneFeatures::GeneFeatures(std::string GTFfilename,std::string id_attribute)
{
    this->total_exon = 0;
    read_features(GTFfilename,id_attribute);
    
}

GeneFeatures::~GeneFeatures(){
    chrom_itvTree_Dict_itr it;
    
    for (it=cds_exon_idx_plus.begin(); it != cds_exon_idx_plus.end(); it++) {
        std::map<int,IntervalTree *> tmp = it->second;
        std::map<int,IntervalTree *>::iterator tmp_itr ;
        //std::cout << tmp.size() << std::endl;
        for (tmp_itr = tmp.begin(); tmp_itr != tmp.end(); tmp_itr ++) {
            delete tmp_itr->second;
        }
        
    }
    for (it=cds_exon_idx_minus.begin(); it != cds_exon_idx_minus.end(); it++) {
        std::map<int,IntervalTree *> tmp = it->second;
        std::map<int,IntervalTree *>::iterator tmp_itr ;
        
        for (tmp_itr = tmp.begin(); tmp_itr != tmp.end(); tmp_itr ++) {
            delete tmp_itr->second;
        }
        
    }
    
}

void GeneFeatures::build_tree(std::map<std::string, std::map<std::string,Gene> > temp_plus, std::map<std::string, std::map<std::string,Gene> > temp_minus)
{
    std::vector<Interval> itemlist;
    std::vector<Interval> sublist;
    std::map<std::string, std::map<std::string,Gene> >::iterator it;
    std::map<std::string,Gene> tmp;
    std::map<std::string,Gene>::iterator tmp_itr;
    gene_exon_Dict_It exon_str_itr;
    
    std::string chr,gid;//,start_ss,end_ss;

    int g_idx =-1;
    int e_idx =0;
    size_t i ;
    int  cur_bin_id,start_bin_id, js, je, k;//, buket_size;

    for (it = temp_plus.begin(); it != temp_plus.end(); it++) {
        chr = it->first;
        tmp = it->second;
        itemlist.clear();
        
        for (tmp_itr = tmp.begin(); tmp_itr != tmp.end(); tmp_itr++) {
            gid = tmp_itr->first;
            features.push_back(gid); //save gene name
            g_idx +=1;
            Gene tmp_gene = tmp_itr->second;
            tmp_gene.get_others();
            
            //int gene_len = (int) tmp_gene.max_stop - tmp_gene.min_start + 1;
            //std::cout << gene_len << std::endl;
            
            std::vector<int> gene_base_pos ;
            
            for (i=0; i< tmp_gene.exons.size(); i++) {
                int st = (int) tmp_gene.exons[i].first - tmp_gene.min_start;
                int end = (int) tmp_gene.exons[i].second - tmp_gene.min_start;
                for (int j=st; j<=end; j++) {
                    gene_base_pos.push_back(j);
                }
            }
            gene_lengths.push_back(tmp_gene.gene_actual_len);
            gene_starts.push_back(tmp_gene.min_start);
            gene_ends.push_back(tmp_gene.max_stop);

            std::vector<int> percentile_list;
            std::sort(gene_base_pos.begin(),gene_base_pos.end());
            
            for (int j=0; j<=100;j++) {
                float kk = (tmp_gene.gene_actual_len - 1) * j/100.0;
                float f = floor(kk);
                float c = ceil(kk);
                if (f == c){
                    percentile_list.push_back( gene_base_pos[kk]);
                }
                else{
                    float d0 = gene_base_pos[int(f)] * (c-kk);
                    float d1 = gene_base_pos[int(c)] * (kk-f);
                    percentile_list.push_back(int(round(d0+d1)));
                }
            }
            gene_percentile_list.push_back(percentile_list);
            
            for (i=0; i < tmp_gene.cds.size(); i++) {
                this->total_exon ++ ;
                //std::cout<< this->total_exon << std::endl;
                itemlist.push_back(Interval(g_idx,e_idx,tmp_gene.cds[i].first,tmp_gene.cds[i].second,CDS));
                e_idx ++;
            }
            for (i=0; i < tmp_gene.utr5.size(); i++) {
                this->total_exon ++ ;
                itemlist.push_back(Interval(g_idx,e_idx,tmp_gene.utr5[i].first,tmp_gene.utr5[i].second,UTR5));
                e_idx ++;
            }

            for (i=0; i < tmp_gene.utr3.size(); i++) {
                this->total_exon ++ ;
                //std::cout << chr << "\t" << tmp_gene.utr3[i].first << "\t" << tmp_gene.utr3[i].second << "\t" << gid << std::endl;
                itemlist.push_back(Interval(g_idx,e_idx,tmp_gene.utr3[i].first,tmp_gene.utr3[i].second,UTR3));
                e_idx ++;
            }
            for (i=0; i < tmp_gene.intron.size(); i++) {
                itemlist.push_back(Interval(g_idx,-1,tmp_gene.intron[i].first,tmp_gene.intron[i].second,INTRON));
            }
            itemlist.push_back(Interval(g_idx,-1,tmp_gene.itg1k[0].first,tmp_gene.itg1k[0].second,ITGUP1K));
            itemlist.push_back(Interval(g_idx,-1,tmp_gene.itg1k[1].first,tmp_gene.itg1k[1].second,ITGDN1K)) ;
            
        }
        
        std::sort(itemlist.begin(),itemlist.end(),itv_comp) ; //key=operator.attrgetter('start'));
        //quick_sort(itemlist,0,itemlist.size()-1);
        start_bin_id = itemlist[0].start/BIN_SIZE;
        js = 0 ;
        je = 0 ;
        k = 0 ;
        
        
        for (i=0; i < itemlist.size(); i++) {
            cur_bin_id = itemlist[i].start/BIN_SIZE;
            //std::cout << cur_bin_id   << std::endl;
            
            if (cur_bin_id == start_bin_id) {
                je += 1;
            }
            else {
                //buket_size = (int)sqrt(je - js) + 1;
                //std::cout << buket_size << std::endl;
                
                sublist = std::vector<Interval>(itemlist.begin()+js,itemlist.begin() + je);

                cds_exon_idx_plus[chr][start_bin_id] = new IntervalTree(sublist);
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
            cds_exon_idx_plus[chr][start_bin_id] = new IntervalTree(sublist);
        //print("tree depth = " + str(cds_exon_idx_plus[chr][start_bin_id].get_depth())+ "\n")
            k+=1;
        }
        
        
    }
    //std::cout << " build minus strand."  << std::endl;
    for (it = temp_minus.begin(); it != temp_minus.end(); it++) {
        //std::cout << "negative strand " << std::endl;
        chr = it->first;
        tmp = it->second;
        itemlist.clear();
        
        for (tmp_itr = tmp.begin(); tmp_itr != tmp.end(); tmp_itr++) {
            gid = tmp_itr->first;
            //std::cout << "GID  " << gid << std::endl;
            Gene tmp_gene = tmp_itr->second;
            tmp_gene.get_others();
            
            features.push_back(gid); //save gene name
            g_idx +=1;
            
            //std::string bs_string = "";
            //int gene_len = (int) tmp_gene.max_stop - tmp_gene.min_start + 1;
            
            //std::vector<std::bitset<100> > bs_list;
            std::vector<int> gene_base_pos ;
            
            /*for (int j=0; j< gene_len; j+=100) {
                std::bitset<100> bs;
                bs_list.push_back(bs);
            }*/
            
            for (i=0; i< tmp_gene.exons.size(); i++) {
                int st = (int) tmp_gene.exons[i].first - tmp_gene.min_start;
                int end = (int) tmp_gene.exons[i].second - tmp_gene.min_start;
                //int size = end - st;
                
                for (int j=st; j<=end; j++) {
                    gene_base_pos.push_back(j);
                }
                
            }
            gene_lengths.push_back(tmp_gene.gene_actual_len);
            gene_starts.push_back(tmp_gene.min_start);
            gene_ends.push_back(tmp_gene.max_stop);

            std::vector<int> percentile_list;
            std::sort(gene_base_pos.begin(),gene_base_pos.end(),reverse_ord_func);
            
            for (int j=0; j<=100;j++) {
                float kk = (tmp_gene.gene_actual_len - 1) * j/100.0;
                float f = floor(kk);
                float c = ceil(kk);
                if (f == c){
                    percentile_list.push_back( gene_base_pos[kk]);
                }
                else{
                    float d0 = gene_base_pos[int(f)] * (c-kk);
                    float d1 = gene_base_pos[int(c)] * (kk-f);
                    percentile_list.push_back(int(round(d0+d1)));
                }
            }
            gene_percentile_list.push_back(percentile_list);
            
            
            for (i=0; i < tmp_gene.cds.size(); i++) {
                total_exon += 1;
                itemlist.push_back(Interval(g_idx,e_idx,tmp_gene.cds[i].first,tmp_gene.cds[i].second,CDS));
                e_idx ++;
            }
            for (i=0; i < tmp_gene.utr5.size(); i++) {
                total_exon += 1;
                itemlist.push_back(Interval(g_idx,e_idx,tmp_gene.utr5[i].first,tmp_gene.utr5[i].second,UTR5));
                e_idx ++;
            }
            for (i=0; i < tmp_gene.utr3.size(); i++) {
                total_exon += 1;
                itemlist.push_back(Interval(g_idx,e_idx,tmp_gene.utr3[i].first,tmp_gene.utr3[i].second,UTR3));
                e_idx ++;
            }
            for (i=0; i < tmp_gene.intron.size(); i++) {
                itemlist.push_back(Interval(g_idx,-1,tmp_gene.intron[i].first,tmp_gene.intron[i].second,INTRON));
            }
            itemlist.push_back(Interval(g_idx,-1,tmp_gene.itg1k[0].first,tmp_gene.itg1k[0].second,ITGUP1K));
            itemlist.push_back(Interval(g_idx,-1,tmp_gene.itg1k[1].first,tmp_gene.itg1k[1].second,ITGDN1K)) ;
            
        }
        
        //std::sort(key=operator.attrgetter('start'));
        //quick_sort(itemlist,0,itemlist.size()-1);
        std::sort(itemlist.begin(),itemlist.end(),itv_comp) ;
        start_bin_id = itemlist[0].start/BIN_SIZE;
        js = 0 ;
        je = 0 ;
        k = 0 ;
        
        for (i=0; i < itemlist.size(); i++) {
            cur_bin_id = itemlist[i].start/BIN_SIZE;
            if (cur_bin_id == start_bin_id) {
                je += 1;
            }
            else {
                //buket_size = (int)sqrt(je - js) + 1;
               // std::cout << buket_size << std::endl;
                
                
                sublist = std::vector<Interval>(itemlist.begin()+js,itemlist.begin() + je);
                //cds_exon_idx_minus[chr][start_bin_id] = new IntervalTree(sublist,16,buket_size,-1,-1,buket_size);
                cds_exon_idx_minus[chr][start_bin_id] = new IntervalTree(sublist);
                k ++;
                start_bin_id = cur_bin_id ;
                js = je ;
                je ++;
            }
        }
        
        
        if (js != je) {
            //buket_size = (int) sqrt(je - js) + 1;
            //std::cout << buket_size << std::endl;
            
            sublist = std::vector<Interval>(itemlist.begin()+js,itemlist.begin() + je);
            //cds_exon_idx_minus[chr][start_bin_id] = new IntervalTree(sublist, 16, buket_size,-1,-1,buket_size );
            cds_exon_idx_minus[chr][start_bin_id] = new IntervalTree(sublist);
            //print("tree depth = " + str(cds_exon_idx_plus[chr][start_bin_id].get_depth())+ "\n")
            k+=1;
        }
        
        
    }

}
//Reading & processing annotation files
void GeneFeatures::read_features(std::string gff_filename, std::string id_attribute)
{

    //dict of dicts since the builtin type doesn't support it for some reason
    std::map<std::string, std::map<std::string, Gene> > temp_plus ;
    std::map<std::string, std::map<std::string, Gene> > temp_minus ;
    std::map<std::string, std::map<std::string, Gene> >::iterator tmp_itr;
    std::map<std::string, Gene>::iterator id_itr;

    bool matched = false;
    //int k = 0;
    int i =  0;
    int counts = 0 ;
    int line_no = 0;
    int start = -1;
    int end = -1;
    std::size_t pos,cur_pos;
    //std::string left_str,sub_str;
    
    std::ifstream input; //(gff_filename);
    
    try{
        input.open (gff_filename, std::ifstream::in);
     
    while(! input.eof()){
    
        std::string line,chrom,source,feature,start_ss,end_ss,score,strand,frame,attributeStr;
        std::stringstream ss;
        std::string id = "";

        if (! std::getline(input,line)){
            break;
        }
        
        line_no ++;
        
        if (line == "\n" || !line.compare(0,1,"#")) {
            continue;
        }
       
        ss << line;
        std::getline(ss,chrom,'\t');
        
        std::getline(ss,source,'\t');
        std::getline(ss,feature,'\t');
        std::getline(ss,start_ss,'\t');
        std::getline(ss,end_ss,'\t');
        std::getline(ss,score,'\t');
        std::getline(ss,strand,'\t');
        std::getline(ss,frame,'\t');
        std::getline(ss,attributeStr,'\t');
        
        //std::cout << strand << std::endl;
        try{
        start = std::stol(start_ss);
        end = std::stol(end_ss);
        }
        catch (const std::invalid_argument& ia) {
            std::cerr << "Invalid argument: " << ia.what() << '\n';
            std::exit(1);

        }
        
        cur_pos = 0;
        //std::cout <<attributeStr << std::endl;
        while (cur_pos < attributeStr.length()) {
            std::size_t next_pos = attributeStr.find(";",cur_pos);
            if (next_pos !=std::string::npos) {
                
                std::string tok = attributeStr.substr(cur_pos,(next_pos - cur_pos));
                //std::cout << tok << std::endl;
                tok = tok.substr(tok.find_first_not_of(' '));
                pos = tok.find('=');
                if (pos == std::string::npos) {
                    pos = tok.find(' ');
                }
                std::string key = tok.substr(0,pos);
                key = key.substr(0,key.find(' '));
                
                std::size_t pos_stop = 1 + pos + tok.substr(pos+1).find_first_not_of(' ');
                std::string val = (tok[pos_stop] == '"') ? tok.substr(pos_stop+1, (tok.length() - (pos_stop+2))) : tok.substr(pos_stop, (tok.length() - (pos_stop+1)));
                
                if (key == id_attribute) {
                    id = val;
                    matched = true;
                    break;
                }
                cur_pos = next_pos + 1;
            }
            else {
                break;
            }
            
        }
        
        if (!matched) {
        
            std::cout << "Failure parsing GFF attribute line." << std::endl;
            exit (EXIT_FAILURE);
        }
        
        if (id == "") {
            continue;
        }
        if (feature == "stop_codon") {
            if (strand == "+" ){
                tmp_itr = temp_plus.find(chrom);
                if (tmp_itr != temp_plus.end()) {
                    id_itr = temp_plus[chrom].find(id);
                    if (id_itr != temp_plus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        
                        (*g).set_stop_codon(start,end);
                    }
                    else{
                        Gene g (id,strand);
                        g.set_stop_codon(start,end);
                        temp_plus[chrom].insert(std::pair<std::string,Gene>(id,g));
                    }
                }
                else{
                    Gene g(id,strand);
                    g.set_stop_codon(start,end);
                    std::map<std::string, Gene>  gene_id_map ;
                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                    temp_plus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    
                }

            }
                    
            if (strand == "-" ) {
                tmp_itr = temp_minus.find(chrom);
                if (tmp_itr != temp_minus.end()) {
                    id_itr = temp_minus[chrom].find(id);
                    if (id_itr != temp_minus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        (*g).set_stop_codon(start,end);
                    }
                    else{
                        Gene g(id,strand);
                        g.set_stop_codon(start,end);
                        std::map<std::string, Gene>  gene_id_map ;
                        temp_minus[chrom].insert(std::pair<std::string,Gene>(id,g));
                        
                    }
                }
                else{
                    Gene g(id,strand);
                    g.set_stop_codon(start,end);
                    std::map<std::string, Gene> gene_id_map ;
                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));

                    temp_minus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    
                    
                }

            }

        }
        if (feature == "CDS" ){
            if (strand == "+" ){
                tmp_itr = temp_plus.find(chrom);
                if (tmp_itr != temp_plus.end()) {
                    id_itr = temp_plus[chrom].find(id);
                    if (id_itr != temp_plus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        //(id_itr->second).add_cds(start,end);
                        (*g).add_cds(start,end);
                    }
                    else{
                        Gene g (id,strand);
                        g.add_cds(start,end);
                        temp_plus[chrom].insert(std::pair<std::string,Gene>(id,g));
                    }
                }
                else{
                    Gene g(id,strand);
                    g.add_cds(start,end);
                    std::map<std::string, Gene>  gene_id_map ;
                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                    temp_plus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    
                }

            }
                    
            if (strand == "-" ) {
                
                tmp_itr = temp_minus.find(chrom);
                if (tmp_itr != temp_minus.end()) {
                    id_itr = temp_minus[chrom].find(id);
                    if (id_itr != temp_minus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        (*g).add_cds(start,end);
                    }
                    else{
                        Gene g(id,strand);
                        g.add_cds(start,end);
                        std::map<std::string, Gene>  gene_id_map ;
                        temp_minus[chrom].insert(std::pair<std::string,Gene>(id,g));
                        
                    }
                }
                else{
                    Gene g(id,strand);
                    g.add_cds(start,end);
                    std::map<std::string, Gene> gene_id_map ;
                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));

                    temp_minus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    
                    
                }

            }
        }
        
        if (feature == "exon" ){
            counts += 1 ;
            if (strand == "+" ){
                tmp_itr = temp_plus.find(chrom);
                if (tmp_itr != temp_plus.end()) {
                    id_itr = temp_plus[chrom].find(id);
                    if (id_itr != temp_plus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        (*g).add_exons(start,end);
                    }
                    else{
                        Gene g (id,strand);
                        g.add_exons(start,end);

                        temp_plus[chrom].insert(std::pair<std::string,Gene>(id,g));
                    }
                }
                else{
                    Gene g (id,strand);
                    g.add_exons(start,end);
                    std::map<std::string, Gene> gene_id_map ;

                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                    temp_plus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    
                }
                
            }
            
            if (strand == "-" ) {

                tmp_itr = temp_minus.find(chrom);
                if (tmp_itr != temp_minus.end()) {
                    
                    id_itr = temp_minus[chrom].find(id);
                    if (id_itr != temp_minus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        (*g).add_exons(start,end);
                    }
                    else{
                        Gene g(id,strand);
                        g.add_exons(start,end);

                        temp_minus[chrom].insert(std::pair<std::string,Gene>(id,g));

                    }
                }
                else{
                    Gene g (id,strand);
                    g.add_exons(start,end);
                    std::map<std::string, Gene> gene_id_map ;
                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));

                    temp_minus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    
                }
            }
                    
        }
    
        i += 1 ;
        //if (i % 100000 == 0 )
        //{
            //sys.stderr.write("%d GTF lines processed.\n" % i);
            //std::cout << i << " GTF lines processed." << std::endl;
        //}
        
    }
    
    
    input.close();
        
    }
    catch(std::ifstream::failure e){
        std::cout << "error in read file " << gff_filename << std::endl;
    }

    if (counts == 0 ){
        std::cout << "Warning: No features of type 'exon' or 'CDS' found in gene GTF file." << std::endl;
    }


    build_tree(temp_plus,temp_minus);

    
}

//find exons of given gene that overlap with the given intervals
//return list of tuples
std::vector<std::pair<int,int> > GeneFeatures::get_exons(std::string chrom,int st,int end,std::string strand)
{
    std::vector<std::pair<int,int> > exons;
    chrom_itvTree_Dict::iterator chrom_it;
    std::map<int,IntervalTree *>::iterator bin_iter;
    std::vector<Interval> fs ;
    std::vector<Interval> temp ;
    size_t i;
    int bin_id ;
    
    //try:
    bin_id = st/BIN_SIZE ;
    
    if (strand == "+" || strand == "."){
        chrom_it = cds_exon_idx_plus.find(chrom);
        if (chrom_it != cds_exon_idx_plus.end()) {
            bin_iter = cds_exon_idx_plus[chrom].find(bin_id);
            if (bin_iter != cds_exon_idx_plus[chrom].end()) {
                fs = (*cds_exon_idx_plus[chrom][bin_id]).find(st,end);
            }
        }
    }
            
    if (strand == "-" || strand == "."){
        chrom_it = cds_exon_idx_minus.find(chrom);
        if (chrom_it != cds_exon_idx_minus.end()) {
            bin_iter = cds_exon_idx_minus[chrom].find(bin_id);
            if (bin_iter != cds_exon_idx_minus[chrom].end()) {
                temp = (*cds_exon_idx_minus[chrom][bin_id]).find(st,end);
                fs.insert(fs.end(),temp.begin(),temp.end());
            }
        }
    }
    
    for(i =0 ; i < fs.size(); i++){
        if (fs[i].type != INTRON && fs[i].type != ITGUP1K && fs[i].type != ITGDN1K){
            int s = fs[i].start;
            int e = fs[i].stop;
            
            if (s < st) {
                s = st;
            }
            if (e > end) {
                e = end;
            }
            exons.push_back(std::make_pair(s,e));
        }
    }

    //std::sort(exons.begin(),exons.end());
    return exons;
}

std::vector<std::string> GeneFeatures::getFeatures() {
    return features ;
}

int GeneFeatures::get_start(int g)
{
    if ((size_t) g < gene_starts.size() && g >=0) {
        return gene_starts[g];
    }
    else{
        return -1;
    }
}

int GeneFeatures::get_stop(int g)
{
    if ((size_t)g < gene_ends.size() && g >=0) {
        return gene_ends[g];
    }
    else{
        return -1;
    }
}

std::string GeneFeatures::get_name(int g)
{
    if ((size_t)g < features.size() && g >=0) {
        return features[g];
    }
    else{
        return "";
    }

}
int GeneFeatures::get_numofgenes(){
    return features.size();
}

int GeneFeatures::exist_in_percentile_list(int gene,int pos){
    int idx = (int)pos - gene_starts[gene];
    for (int i=0; i<=100; i++) {

        if (idx == gene_percentile_list[gene][i]) {
            return i;
        }

    }
    return -1;
}



std::map<int,int> GeneFeatures::Gene_annotation(std::string chrom, std::vector<std::pair<int,int> > itv_list,std::string strand,std::vector<int> * mapped_exons )
{
    std::vector<Interval> fs ;
    chrom_itvTree_Dict::iterator chrom_it;
    std::map<int,IntervalTree *>::iterator bin_iter;
    
	size_t i,j;
    int bin_id_s, bin_id_e;
    
    for (i=0; i < itv_list.size(); i ++) {
        bin_id_s = itv_list[i].first/BIN_SIZE;
        bin_id_e = itv_list[i].second/BIN_SIZE;
        
        if (strand == "+" || strand == ".") {
            chrom_it = cds_exon_idx_plus.find(chrom);
            if (chrom_it !=  cds_exon_idx_plus.end()) {
                
                bin_iter = cds_exon_idx_plus[chrom].find(bin_id_s);
                if (bin_iter != cds_exon_idx_plus[chrom].end()  )
                {
                    fs = (*cds_exon_idx_plus[chrom][bin_id_s]).find(itv_list[i].first,itv_list[i].second);
                    if (bin_id_s != bin_id_e) {
                        bin_iter = cds_exon_idx_plus[chrom].find(bin_id_e);
                        if (bin_iter != cds_exon_idx_plus[chrom].end()  ){
                        std::vector<Interval>  tmp = (*cds_exon_idx_plus[chrom][bin_id_e]).find(itv_list[i].first,itv_list[i].second);
                        fs.insert(fs.end(),tmp.begin(),tmp.end());
                        }
                    }
                }
            }
        }

        if (strand == "-" or strand == ".") {
            chrom_it = cds_exon_idx_minus.find(chrom);
            if (chrom_it !=  cds_exon_idx_minus.end()) {
                bin_iter = cds_exon_idx_minus[chrom].find(bin_id_s);
                if (bin_iter != cds_exon_idx_minus[chrom].end()  )
                {
                    std::vector<Interval>  tmp = (*cds_exon_idx_minus[chrom][bin_id_s]).find(itv_list[i].first,itv_list[i].second);
                    fs.insert(fs.end(),tmp.begin(),tmp.end());

                    if (bin_id_s != bin_id_e) {
                        bin_iter = cds_exon_idx_minus[chrom].find(bin_id_e);
                        if (bin_iter != cds_exon_idx_minus[chrom].end()  ){
                        std::vector<Interval>  tmp = (*cds_exon_idx_minus[chrom][bin_id_e]).find(itv_list[i].first,itv_list[i].second);
                        fs.insert(fs.end(),tmp.begin(),tmp.end());
                        }
                    }
                }
           }
        }
    }
    
    
    std::map<int,int> gene_type_map;
    std::map<int,int> gene_ovp_len_map;
    
    std::map<int,int> None_gene_type_map;
    std::map<int,int>::iterator gene_type_map_itr;
    int min_type = 6;
    for (i=0; i < fs.size(); i++) {
        if (fs[i].type <= UTR3) {
            
                int ovp_len = 0;
                for(j=0; j < itv_list.size(); j ++) {
                    if (itv_list[j].first < fs[i].stop && itv_list[j].second > fs[i].start) {
                        int ovp_st = itv_list[j].first > fs[i].start ? itv_list[j].first : fs[i].start;
                        int ovp_end = itv_list[j].second > fs[i].stop ? fs[i].stop : itv_list[j].second;
                        ovp_len += ovp_end - ovp_st+1;
                    }
                }
            if (gene_ovp_len_map.find(fs[i].gene) != gene_ovp_len_map.end()) {
                gene_ovp_len_map[fs[i].gene] += ovp_len;
                if (fs[i].type <= min_type) {
                    gene_type_map[fs[i].gene] = fs[i].type;
                    min_type = fs[i].type;
                }
                
            } else {
                gene_ovp_len_map.insert(std::pair<int,int> (fs[i].gene,ovp_len));
                gene_type_map.insert(std::pair<int,int>(fs[i].gene,fs[i].type));
                min_type = fs[i].type;
            }
        }
        else { //not CDS or UTR
            if (None_gene_type_map.find(fs[i].gene) != None_gene_type_map.end())
            {
               if(fs[i].type < None_gene_type_map[fs[i].gene])
               {
                   None_gene_type_map[fs[i].gene] = fs[i].type;
               }
            }
            else {
                None_gene_type_map.insert(std::pair<int, int> (fs[i].gene,fs[i].type));
            }
        }

         //genes.push_back(std::pair<int,int>(fs[i].type,fs[i].gene));
        if (fs[i].exon != -1) {
            mapped_exons->push_back(fs[i].exon);
        }
    }
    
    if (gene_type_map.size() == 0) {
        return None_gene_type_map;
    }
    else {
        std::map<int,int> res_gene_type_map;
        int max_ovp_len = 0;
        for (auto& kv : gene_ovp_len_map) {
            if( max_ovp_len < kv.second) { max_ovp_len = kv.second; }
        }
        for (auto& kv : gene_ovp_len_map) {
            if (gene_ovp_len_map[kv.first] == max_ovp_len) {
                res_gene_type_map.insert(std::pair<int,int> (kv.first,gene_type_map[kv.first]));
            }
        }
        return res_gene_type_map;
    }
    
    return gene_type_map;

}
/*

int main() {
    std::string filename = "test.gtf";
    std::string id_attr = "gene_id";
   
    std::vector<std::pair<int,int> > itv_list ;
    chr_ITV exp ;
    exp.chrom = "chr4";
    int start = 1048489;
    int end = 1049900;
    std::pair<int,int> p (start,end);
    
    itv_list.push_back(p);
    
    std::cout << "start to build tree " << std::endl;
    GeneFeatures gIdx (filename,id_attr);
    std::cout << "after  build tree " << std::endl;
    
    std::cout << "total exon " << gIdx.total_exon << std::endl;
    std::cout << "total gene " << gIdx.features.size() << std::endl;

   // for (int i=0; i < itv_list.size(); i++) {
    //    std::cout << itv_list[i].start << std::endl;
    //}
    std::vector<int> *exons = new std::vector<int>();
    std::vector<int> res = gIdx.Gene_annotation("chr4",itv_list,".",exons);
    std::vector<std::pair<int,int> > exon_list = gIdx.get_exons("chr4",start,end,".");
    for (auto& p : exon_list){
        std::cout << p.first << "\t" << p.second << std::endl;
    }
    std::vector<std::vector<int> > gene_percentile_list = gIdx.gene_percentile_list;
    
    for (int i=0; i<gene_percentile_list.size(); i++) {
        std::cout << i << std::endl;
        std::vector<int> perc_list = gene_percentile_list[i];
        for (int j=0; j<perc_list.size(); j++) {
            std::cout << perc_list[j] << "\t";
        }
        std::cout << "\n";
    }
    
    std::cout << exons->size() << std::endl;
    
    for (int i=0;i<res.size(); i++) {
        std::cout << res[i] << std::endl;
        //std::cout << exons->operator[](i) << std::endl;
        
    }
    delete exons;
    

    bool test_bool = false;
    std::cout << test_bool << std::endl;
    
}*/
