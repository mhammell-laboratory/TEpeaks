//
//  EM_TEpeaks.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/22/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____EM_TEpeaks__
#define __TEToolkit_c____EM_TEpeaks__

#include <stdio.h>
#include <string>
#include "Parser.h"
#include "ShortRead.h"


extern "C" {

    
    int run_EM_TEpeaks(opt_t options,ShortRead * treat, ShortRead * control,std::string fname);
}


#endif /* defined(__TEToolkit_c____EM_TEpeaks__) */
