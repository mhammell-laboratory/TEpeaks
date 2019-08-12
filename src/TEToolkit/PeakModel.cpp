//
//  PeakModel.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 4/29/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//
// This class is modified from MACS2.


//#include <map>
#include "PeakModel.h"
#include <math.h>
#include <limits.h>
//#include <vector>
#include <numeric>
#include <algorithm>
#include "boost/foreach.hpp"

#include "zeroin.h"
#include "myLog.h"

#include "Constants.h"



PeakModel::PeakModel(ShortRead *treatment, opt_t options, int max_pairnum )
{
    this->treatment = treatment;
    
    gz = options.gsize;
    umfold = options.umfold;
    lmfold = options.lmfold;
    tag_expansion_size = 10;
    this->NotEnoughPairs = false;
    
    //opt.tsize| test 10bps. The reason is that we want the best 'lag' between left & right cutting sides. A tag will be expanded to 10bps centered at cutting point.
    bw = options.bandwidth;
    //info  = opt.info;
    //debug = opt.debug;
    //warn  = opt.warn;
    //error = opt.warn;
    
    this->max_pairnum = max_pairnum;

    if(!build())
    {
        this->NotEnoughPairs = true;
    }
    
}

PeakModel::~PeakModel(){}

bool PeakModel::build ()
{
    //Build the model.

    //prepare d, scan_window, plus_line,
    //minus_line and shifted_line to use.

    //int num_paired_peakpos, num_paired_peakpos_remained, num_paired_peakpos_picked;
    
    peaksize = 2 * bw;
    debug("peak size = " + std::to_string(peaksize));
    min_tags = int(std::round(1.0 * treatment->total_uniq * lmfold * peaksize / gz /2));
    // mininum unique hits on single strand
    max_tags = int(std::round(1.0 * treatment->total_uniq * umfold * peaksize / gz /2));
    // maximum unique hits on single strand

    // use treatment data to build model
    info("#3 looking for paired plus/minus strand peaks...");
    PairedPeaks_Dict paired_peakpos;
    __paired_peaks (paired_peakpos);
    
    // select up to 1000 pairs of peaks to build model
    int num_paired_peakpos = 0;
//ADD    int num_paired_peakpos_remained = max_pairnum;
    int num_paired_peakpos_picked = 0;
    //select only num_paired_peakpos_remained pairs.
    
    BOOST_FOREACH(PairedPeaks_Dict::value_type pp, paired_peakpos)
    {
        num_paired_peakpos += pp.second.size() ;
    }
    
        // TL: Now I want to use everything

    num_paired_peakpos_picked = num_paired_peakpos;

    info("#3 number of paired peaks: " + std::to_string(num_paired_peakpos));
    if (num_paired_peakpos < 100){
        char buffer[500];
        sprintf(buffer, "Too few paired peaks (%d) so I can not build the model! Broader your MFOLD range parameter may erase this error. If it still can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed number instead.",num_paired_peakpos);
        
        error(std::string(buffer));
        error("#3 Process for pairing-model is terminated!");
        return false;
    }
    else {
        if(num_paired_peakpos < max_pairnum) {
            char buffer[500];
            sprintf(buffer, "Fewer paired peaks (%d) than %d! Model may not be build well! Lower your MFOLD parameter may erase this warning. Now I will use %d pairs to build model!",num_paired_peakpos,max_pairnum,num_paired_peakpos_picked);
            warn(std::string(buffer));
        }
    }
    char buffer[200] ;
    sprintf(buffer,"Use %d pairs to build the model.",num_paired_peakpos_picked );
    //debug(std::string(buffer));
    
    __paired_peak_model(paired_peakpos);
    
    return true;
        
}


void PeakModel::__paired_peak_model ( PairedPeaks_Dict paired_peakpos)
{
        /*Use paired peak positions and treatment tag positions to build the model.

        Modify (d, model_shift size and scan_window size. and extra, plus_line, minus_line and shifted_line for plotting).
        */

    //std:vector<int> paired_peakpos_chrom;
    
    std::vector<double> plus_data, minus_data;
        
    int window_size = 1 + 2 * peaksize + tag_expansion_size;
    
    std::vector<int> plus_start(window_size,0);     // for fast pileup
    std::vector<int> plus_end(window_size,0);      // for fast pileup
    std::vector<int> minus_start(window_size,0);    // for fast pileup
    std::vector<int> minus_end(window_size,0);       // for fast pileup
    //std::vector<int> plus_line(window_size,0);
    //std::vector<int> minus_line(window_size,0);
    for (int k=0; k < window_size; k++) {
        plus_line[k] = 0;
        minus_line[k] = 0;
    }
    
    //info("start model_add_line...");
    //std::vecotr<int> chroms = paired_peakpos.keys();
        
    //BOOST_FOREACH( PairedPeaks_Dict::value_type pp, paired_peakpos)
    for (auto pp : paired_peakpos)
    {
        //std:vector<int>  paired_peakpos_chrom = paired_peakpos[pp.first];
        
        std::vector<int> tags_plus = treatment->get_locations_by_chr_uniq( pp.first, 1);
        std::vector<int> tags_minus  = treatment->get_locations_by_chr_uniq( pp.first, 0 );
        
        // every paired peak has plus line and minus line

        __model_add_line (pp.second, tags_plus, plus_start, plus_end);
        __model_add_line (pp.second, tags_minus, minus_start, minus_end) ;

    }
    __count ( plus_start, plus_end, plus_line );
    __count ( minus_start, minus_end, minus_line );

    info("start X-correlation...");
    // Now I use cross-correlation to find the best d
        
    // normalize first
    double minus_avg = (std::accumulate(minus_line.begin(),minus_line.begin() + window_size,0))/window_size;
    double plus_avg = (std::accumulate(plus_line.begin(),plus_line.begin() + window_size,0))/window_size;
    
    double minus_std = standard_deviation(minus_line,minus_avg);
    double plus_std = standard_deviation(plus_line,plus_avg);
    
    for (int k=0; k < window_size; k++) {
        minus_data[k] = (minus_line[k] - minus_avg)/(minus_std * window_size);
        plus_data[k] = (plus_line[k] - plus_avg)/(plus_std * window_size);
    }

    //cross-correlation
    this->ycorr = correlate(minus_data,plus_data);
    std::vector<double> ycorr_2peaksize ;
    for(int k= window_size-peaksize;k < window_size + peaksize; k++)
    {
        ycorr_2peaksize.push_back(ycorr[k]);
    }
    
    
    linspace(int(ycorr_2peaksize.size())/2 *-1, int(ycorr_2peaksize.size())/2, int(ycorr_2peaksize.size()), this->xcorr);

    //smooth correlation values to get rid of local maximums from small fluctuations.
    smooth(ycorr_2peaksize, 11,"flat"); // window size is by default 11.

    //all local maximums could be alternative ds.
    //int i_l_max[ycorr_2peaksize.size()];
    
    //i_l_max[0] = 0;
    //i_l_max[ycorr_2peaksize.size()-1] = 0;
    
    //std::vector<double> tmp_cor_alternative_d;
    //std::vector<double> tmp_alternative_d;
    std::vector<double> cor_alternative_d;
    std::vector<int> alternative_d;
    int max_corr = INT_MIN;
    
    
    for (size_t k= 1; k < ycorr_2peaksize.size() -1 ; k++) {
        
            if ((ycorr_2peaksize[k] > ycorr_2peaksize[k-1]) && (ycorr_2peaksize[k] > ycorr_2peaksize[k+1])) {
                //i_l_max[k] = 1;
                
                if (xcorr[k] >0) {
                    cor_alternative_d.push_back(ycorr_2peaksize[k]);
                    alternative_d.push_back(int(xcorr[k]));
                    
                    if (ycorr_2peaksize[k] > max_corr) {
                        max_corr = ycorr_2peaksize[k];
                        this->d = xcorr[k];
                    }
                }
            }
        this->ycorr[k] = ycorr_2peaksize[k];
        
    }
    
    this->ycorr[0] = ycorr_2peaksize[0];
    //i_l_max = np.r_[False, ycorr[1:] > ycorr[:-1]] & np.r_[ycorr[:-1] > ycorr[1:], False];
    //tmp_cor_alternative_d = ycorr[ i_l_max ];
   // tmp_alternative_d = xcorr[ i_l_max ];
    //cor_alternative_d =  tmp_cor_alternative_d [ tmp_alternative_d > 0 ];
   // alternative_d = map( int, tmp_alternative_d[ tmp_alternative_d > 0 ] )
        
        // best cross-correlation point
    //d = xcorr[ np.where( ycorr== max( cor_alternative_d ) )[0][0] ];
        //d = xcorr[np.where(ycorr==max(ycorr))[0][0]] #+tag_expansion_size

        // get rid of the last local maximum if it's at the right end of curve.
        
    if(cor_alternative_d.size() > 0 )
    {
        warn("No proper d can be found! Tweak --mfold?");
    }
        
    scan_window = std::max(d,tag_expansion_size)*2;
    info("end of X-cor");
        
    //return true;
}
                            
void PeakModel::__model_add_line ( std::vector<int> pos1, std::vector<int> tags, std::vector<int>& start, std::vector<int>& end)
{
        /*"""Project each pos in pos2 which is included in  [pos1-peaksize,pos1+peaksize] to the line.

        pos1: paired centers -- array.array
        pos2: tags of certain strand -- a numpy.array object
        line: numpy array object where we pileup tags

        */
        
    int i1 = 0;                  // index for pos1
    int i2 = 0;                  // index for pos2
    int i2_prev = 0;             //index for pos2 in previous pos1
                                // [pos1-peaksize,pos1+peaksize]
                                // region
    int i1_max = pos1.size();
    int i2_max = tags.size();
    
    bool flag_find_overlap = false;
    int max_index = start.size() - 1;
    
    int psize_adjusted1 = peaksize + tag_expansion_size / 2; //# half window

    while( i1 < i1_max and i2< i2_max){
        
            int p1 = pos1[i1];
            int p2 = tags[i2];
                
            if (p1 - psize_adjusted1 > p2){ // move pos2
                i2 += 1;
            }
            else {
                
            if (p1 + psize_adjusted1 < p2){ // move pos1
                i1 += 1;
                i2 = i2_prev;    // search minus peaks from previous index
                flag_find_overlap = false;
            }
            else{              // overlap!
                if( not flag_find_overlap)
                {
                    flag_find_overlap = true;
                    i2_prev = i2; // only the first index is recorded
                }
                int s = std::max(p2 - tag_expansion_size/2 - p1 + psize_adjusted1, 0);
                start[s] += 1;
                int e = std::min(p2 + tag_expansion_size/2 - p1 + psize_adjusted1, max_index);
                end[e] -= 1;

                i2+=1;
                }
            }
    }
    
}
                        
void PeakModel::__count ( std::vector<int> start, std::vector<int> end, std::vector<int>& line )
{
    int pileup = 0;
    
    for(size_t i=0;i< line.size(); i++)
    {
        pileup += start[i] + end[i];
        line[i] = pileup;
    }
}


void PeakModel::__paired_peaks (PairedPeaks_Dict &paired_peaks_pos)
{
    /* Call paired peaks from fwtrackI object. Return paired peaks center positions.*/
    
    //int i;
    std::vector<std::string> chrs;
    //std::string chrom;
   // PairedPeaks_Dict paired_peaks_pos;
    
    std::vector<int> plus_tags;
    std::vector<int> minus_tags;
    
    for(auto v : treatment->chromNameToId){
        
        //debug("Chromosome: " + v.first );
        
        plus_tags = treatment->get_locations_by_chr_uniq( v.second, 1);
        minus_tags  = treatment->get_locations_by_chr_uniq( v.second, 0 );
        
        debug("before naive find peaks " );
        debug("size = " + std::to_string(treatment->numOfUniqTags_by_chr(v.second,0)));
        
        std::vector<std::pair<int, int> > plus_peaksinfo = __naive_find_peaks ( plus_tags, treatment->numOfUniqTags_by_chr(v.second,0), 1 );
        
        debug("Number of unique tags on + strand: " + std::to_string( treatment->numOfUniqTags_by_chr(v.second,0)) );
        debug("Number of peaks in + strand: " + std::to_string( plus_peaksinfo.size() ) );
        
        std::vector<std::pair<int, int> > minus_peaksinfo = __naive_find_peaks ( minus_tags, treatment->numOfUniqTags_by_chr(v.second,1), 0 );
        
        //debug("Number of unique tags on - strand: " + std::to_string( treatment->numOfUniqTags_by_chr(v.second,1) ) );
        //debug("Number of peaks in - strand: " + std::to_string( minus_peaksinfo.size() ) );
        
        if (plus_peaksinfo.size()==0 || minus_peaksinfo.size()==0){
            debug("Chromosome " + v.first + " is discarded!" );
            continue;
        }
        else{
            paired_peaks_pos.insert(std::pair<int, std::vector<int> > (v.second, __find_pair_center (plus_peaksinfo, minus_peaksinfo) ) );
            //debug("Number of paired peaks: " + std::to_string(paired_peaks_pos[v.second].size()));
            
        }
    }
    //return paired_peaks_pos;
}

std::vector<int> PeakModel::__find_pair_center (std::vector<std::pair<int,int> > pluspeaks, std::vector<std::pair<int,int> > minuspeaks)
{
    int ip = 0;                  // index for plus peaks
    int im = 0;                  // index for minus peaks
    int im_prev = 0;             // index for minus peaks in previous plus peak
    std::vector<int> pair_centers ;
            
    int ip_max = pluspeaks.size();
    int im_max = minuspeaks.size();
            
    bool flag_find_overlap = false;
    
    while (ip < ip_max and im < im_max){
        int pp = pluspeaks[ip].first;// # for (peakposition, tagnumber in peak)
        int pn = pluspeaks[ip].second;
        int mp = minuspeaks[im].first;
        int mn = minuspeaks[im].second;
        
        if (pp-peaksize > mp){ //# move minus
            im += 1;
        }
        else {
            if( pp + peaksize < mp){ // # move plus
                ip += 1;
                im = im_prev;    // search minus peaks from previous index
                flag_find_overlap = false;
            }
            else{              // overlap!
                if(not flag_find_overlap){
                    flag_find_overlap = true;
                    im_prev = im;// # only the first index is recorded
                }
                if(1.0 * pn/mn < 2 and 1.0 * pn/mn > 0.5) // # number tags in plus and minus peak region are comparable...
                {
                    if (pp < mp){
                        pair_centers.push_back((pp+mp)/2);
                    }
                        //debug ( "distance: %d, minus: %d, plus: %d" % (mp-pp,mp,pp))
                }
                im += 1;
            }
        }
    }
    return pair_centers;
}

std::vector<std::pair<int, int> > PeakModel::__naive_find_peaks (std::vector<int> taglist, int size, int plus_strand )
{
        /*"""Naively call peaks based on tags counting. if plus_strand == 0, call peak on minus strand. Return peak positions and the tag number in peak region by a tuple list [(pos,num)].
        """*/
    //long i;
    int pos;
    std::vector<std::pair<int,int> > peak_info;

    std::vector<int> current_tag_list;

    if (size < 2)
    {
        return peak_info;
    }
    
    pos = taglist[0];
    
    current_tag_list.push_back(pos);
    
    debug("list size = " + std::to_string(size));
    for (int i=1; i < size; i++)
    {
        pos = taglist[i];

        if (( pos - current_tag_list[0] + 1 ) > peaksize)
        {
            // call peak in current_tag_list when the region is long enough
            // a peak will be called if tag number is ge min tags.
            debug("pos = " + std::to_string(pos));
            debug("current_tag_list[0] = " + std::to_string(current_tag_list[0]));
            debug("min_tags = " + std::to_string(min_tags));
            debug("max_tags = " + std::to_string(max_tags));
            debug("current tab list size = " + std::to_string(current_tag_list.size()));
            if ((int)current_tag_list.size() >= min_tags && (int)current_tag_list.size() <= max_tags)
            {
                debug("current i = " + std::to_string(i));
                peak_info.push_back( std::pair<int,int>( __naive_peak_pos(current_tag_list,plus_strand), current_tag_list.size()) );
                debug("after naive peak pos");
            }

            current_tag_list.clear();
        }
        current_tag_list.push_back( pos );
    }
    return peak_info;
}

int PeakModel::__naive_peak_pos (std::vector<int> pos_list, int plus_strand )
{
        /*"""Naively calculate the position of peak. plus_strand: 1, plus; 0, minus return the highest peak summit position.
        """*/

  //  int  s, e;
    
    std::vector<int> ss, es;

    int peak_length = pos_list.back() + 1 - pos_list.front() + tag_expansion_size;
    
    debug("peak_length = " + std::to_string(peak_length));

    int start = pos_list[0] - tag_expansion_size/2; // leftmost position of project line

    std::vector<int> horizon_line (peak_length,0);  // the line for tags to be projected

    debug("start in naive peak pos 0");
    
    BOOST_FOREACH(int pos, pos_list)
    {
        ss.push_back( std::max(pos - start - tag_expansion_size/2,0) );
        es.push_back( std::min(pos - start + tag_expansion_size/2,peak_length) );
    }

    std::sort(ss.begin(),ss.end());
    std::sort(es.begin(),es.end());
    
    int pileup = 0;

    int ls =  ss.size();
    int le = es.size();

    int i_s = 0;
    int i_e = 0;
    
    int pre_p = std::min( ss[ 0 ], es[ 0 ] );
    int i=0;
    
    if (pre_p != 0){
        for (;i< pre_p;i++ ){
            horizon_line[ i ] = 0;
        }
    }
    debug("start in naive peak pos 1");
    while( i_s < ls and i_e < le){
        if( ss[ i_s ] < es[ i_e ]){
            int p = ss[ i_s ];
            if (p != pre_p){
                for( i = pre_p;i < p; i++ ){
                    horizon_line[ i ] = pileup;
                }
                pre_p = p;
            }
            pileup += 1;
            i_s += 1;
        }
        else {
                if(ss[ i_s ] > es[ i_e ])
                {
                    int p = es[ i_e ];
                    if (p != pre_p){
                        for (i=pre_p;i< p;i++){
                            horizon_line[ i ] = pileup;
                        }
                        pre_p = p;
                    }
                    pileup -= 1;
                    i_e += 1;
                }
                else{
                        i_s += 1;
                        i_e += 1;
                }
            
            }
    }
    if ( i_e < ls ){
        for(int j=i_e;j<ls;j++){
            int p = es[ i_e ];
            if (p != pre_p){
                for (i=pre_p;i<p;i++){
                    horizon_line[ i ] = pileup;
                }
                pre_p = p;
            }
            pileup -= 1;
        }
    }
    debug("naive peak pos 2");
    std::vector<int> top_pos ;            // to record the top positions. Maybe > 1
    int top_p_num = 0 ;          // the maximum number of projected points
    
    for (int pp = 0;pp< peak_length;pp++){ // find the peak posistion as the highest point
        if (horizon_line[pp] > top_p_num){
            top_p_num = horizon_line[pp];
            top_pos.clear();
            top_pos.push_back(pp);
        }
        else {
            if(horizon_line[pp] == top_p_num){
                top_pos.push_back(pp);
            }
        }
    }
      debug("naive peak pos 3");
    return top_pos[int(top_pos.size()/2)]+start;
}

