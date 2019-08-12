//
//  EMestimate_reads.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/23/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#include <vector>
#include <limits>
#include <iostream>
#include <cmath>

#include "EMestimate_reads.h"
#include "myLog.h"


void EMUpdate(std::vector<double> & means,std::vector<int> multiAlgn_To_multiRead,std::vector<int>multiAlgn_To_peakID, std::vector<double> & peak_reads_Prime,std::vector<int> effLengths,std::vector<double> & meansPrime )
{
    //update peak reads according to means (weight of each multi-alignment)
    int multiAlgn_first_idx = 0;
    //std::vector<double> tmp_peak_reads(peak_reads.size(),0.0);
    
    for (size_t i=0; i< multiAlgn_To_multiRead.size(); i++) {
        //loop through all multi-alignments of each multi-read
        for (int j= multiAlgn_first_idx; j < multiAlgn_To_multiRead[i]; j++) {
            int pid = multiAlgn_To_peakID[j];
            peak_reads_Prime[pid] += means[j];
        }
        multiAlgn_first_idx = multiAlgn_To_multiRead[i];
    }
    //for (int k =0; k < peak_reads.size(); k++) {
    //    std::cout << peak_reads[k] << std::endl;
    //}
    //compute means based on updated peak reads
    multiAlgn_first_idx = 0;
    for (size_t i=0; i < multiAlgn_To_multiRead.size(); i++) {
        double total = 0.0;
        
        //loop through all multi-alignments of each multi-read to get total counts
        for (int j= multiAlgn_first_idx; j < multiAlgn_To_multiRead[i]; j++) {
            int pid = multiAlgn_To_peakID[j];
            total += 1.0 * peak_reads_Prime[pid]/effLengths[pid];
        }
        //
        for (int j= multiAlgn_first_idx; j < multiAlgn_To_multiRead[i]; j++) {
            
            int pid = multiAlgn_To_peakID[j];
            meansPrime[j] = 1.0 *peak_reads_Prime[pid]/(effLengths[pid]*total);
            
        }
        multiAlgn_first_idx = multiAlgn_To_multiRead[i];
    }
    
}

/*
 *compute weight of each multi-alignment using EM algorithm
 *multi-alignments of one multi-read have consecutive ids.
 *multiAlgn_To_multiRead : multi-read -> last index of its multi-alignments
 *
 * means: weight of multi_alignment
 */

void EMestimate_read(std::vector<int> multiAlgn_To_multiRead,std::vector<int> multiAlgn_To_peakID,std::vector<double> & means,std::vector<double> peak_reads, std::vector<double> & peak_reads_Prime, int numItr,std::vector<int> effLengths)
{
    //minimum weight
    //double alpha_limit = 1e-7;
    //minimum changes
    double alpha_diff_cutoff = 0.01;
    double alpha_check_cutoff = 0.01;
    bool converged = false;
    
    //std::vector<double> means(multiAlgn_weight.size(),0.0);
    std::vector<double> meansPrime(means.size(),0.0);
    //std::vector<double> peak_readsPrime(peak_reads.size(),0.0);

    int cur_iter = 0;
    
    while (cur_iter < numItr && !converged)
    {
        cur_iter += 1;
        for (size_t i=0; i < peak_reads.size(); i++) {
            peak_reads_Prime[i] = peak_reads[i];
        }

       // debug("interation " + std::to_string(cur_iter));
    EMUpdate(means,multiAlgn_To_multiRead,multiAlgn_To_peakID,peak_reads_Prime,effLengths,meansPrime);
        
        converged = true;
        double maxRelDiff = -std::numeric_limits<double>::max();
        
        for (size_t i = 0; i < means.size(); ++i) {
            if (meansPrime[i] > alpha_check_cutoff) {
                double relDiff = std::fabs(means[i] - meansPrime[i]) / meansPrime[i];
                
                maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
                if (relDiff > alpha_diff_cutoff) {
                    converged = false;
                }
            }
            means[i] = meansPrime[i];
            meansPrime[i] = 0.0;
        }
        cur_iter ++;
        
        if (converged) {
            //std::cout << "Converged in " << cur_iter << " iterations!" << std::endl;
            info("Converged in " + std::to_string(cur_iter) + " iterations!");
        }
        
    }
    
    //output final reads for each peak including both unique reads and multi-reads
    
    for (size_t i=0; i < peak_reads.size(); i++) {
        peak_reads_Prime[i] = peak_reads[i];
    }
    //update peak reads according to means (weight of each multi-alignment)
    int multiAlgn_first_idx = 0;
    //std::vector<double> tmp_peak_reads(peak_reads.size(),0.0);
    
    for (size_t i=0; i< multiAlgn_To_multiRead.size(); i++) {
        //loop through all multi-alignments of each multi-read
        for (int j= multiAlgn_first_idx; j < multiAlgn_To_multiRead[i]; j++) {
            int pid = multiAlgn_To_peakID[j];
            peak_reads_Prime[pid] += means[j];
        }
        multiAlgn_first_idx = multiAlgn_To_multiRead[i];
    }

    
    
   /* for (int i=0; i < means.size(); i++) {
        if (means[i] < alpha_limit) {
            means[i] = 0.0;
        }
    }*/
 
}

/*
int main() {
 
    std::vector<int> multiAlgn_To_multiRead;
    std::vector<int> multiAlgn_To_peakID;
    std::vector<double> means;
    std::vector<double> peak_reads;
    int numItr = 100;
    std::vector<int> effLengths;
    std::vector<double> peak_reads_Prime(4,0.0);
     
    peak_reads.push_back(4.0);
    peak_reads.push_back(3.0);
    peak_reads.push_back(3.0);
    peak_reads.push_back(3.0);
    
    effLengths.push_back(2);
    effLengths.push_back(2);
    effLengths.push_back(2);
    effLengths.push_back(2);
    
    means.push_back(0.4);
    means.push_back(0.3);
    means.push_back(0.3);
    means.push_back(0.33);
    means.push_back(0.33);
    means.push_back(0.33);
    
    multiAlgn_To_peakID.push_back(0);
    multiAlgn_To_peakID.push_back(1);
    multiAlgn_To_peakID.push_back(2);
    multiAlgn_To_peakID.push_back(1);
    multiAlgn_To_peakID.push_back(2);
    multiAlgn_To_peakID.push_back(3);
    
    multiAlgn_To_multiRead.push_back(3);
    multiAlgn_To_multiRead.push_back(6);
 
    EMestimate_read(multiAlgn_To_multiRead,multiAlgn_To_peakID,means,peak_reads,peak_reads_Prime,numItr,effLengths);
 
    for (int i=0; i < peak_reads_Prime.size(); i++) {
        std::cout << peak_reads_Prime[i]  << std::endl;
    }
    
    for (int i=0; i < means.size(); i++) {
        std::cout << means[i]  << std::endl;
    }
    
 
 
 }
*/
