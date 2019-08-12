//
//  EMestimate_reads.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/23/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____EMestimate_reads__
#define __TEToolkit_c____EMestimate_reads__



#endif /* defined(__TEToolkit_c____EMestimate_reads__) */

extern "C" {

    
    void EMestimate_read(std::vector<int> multiAlgn_To_multiRead,std::vector<int> multiAlgn_To_peakID,std::vector<double> & means,std::vector<double> peak_reads, std::vector<double> & peak_reads_Prime, int numItr,std::vector<int> effLengths);
    
}