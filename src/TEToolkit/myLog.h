//
//  myLog.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 4/29/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____myLog__
#define __TEToolkit_c____myLog__

#include <stdio.h>
#include <string>
// initialize thy logs..
void info(std::string message);
void warn(std::string message);
void debug(std::string message);
void error(std::string message);

class Test{
public:
    int t;
    Test();
    ~Test();
};
#endif /* defined(__TEToolkit_c____myLog__) */
