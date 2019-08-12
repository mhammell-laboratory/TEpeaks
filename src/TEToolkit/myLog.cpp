//
//  myLog.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 4/29/16.
//  Copyright (c) 2016 Ying Jin. All rights reserved.
//



#include "myLog.h"

#include <iostream>
//#include <iomanip>
#include <string>
#include <stdlib.h>
#include <vector>
#include <map>
#include <ctime>


//#include "boost/date_time/posix_time/posix_time.hpp"
//#include <boost/algorithm/string/join.hpp>
//#include "boost/foreach.hpp"
//#include <boost/accumulators/accumulators.hpp>
//#include <boost/accumulators/statistics.hpp>
//#include <boost/range/adaptor/transformed.hpp>



void info(std::string m) {
    // formatting    : time [idx] message \n
    // destinations  : console, file "out.txt" and debug window
    //namespace pt = boost::posix_time;
    std::time_t t = std::time(nullptr);
    
   // pt::ptime now = pt::second_clock::local_time();
    
    //std::string now_str = to_simple_string(now);
    
    std::cout << std::asctime(std::localtime(&t)) << "\tINFO\t" << m << std::endl;
    
}

void warn(std::string m) {
    // formatting    : time [idx] message \n
    // destinations  : console, file "out.txt" and debug window
    std::time_t t = std::time(nullptr);
   // namespace pt = boost::posix_time;
   // pt::ptime now = pt::second_clock::local_time();
    
   // std::string now_str = to_simple_string(now);
    
    std::cout << std::asctime(std::localtime(&t)) << "\tWARNING\t" << m << std::endl;
    
}

void debug(std::string m) {
    // formatting    : time [idx] message \n
    // destinations  : console, file "out.txt" and debug window
    std::time_t t = std::time(nullptr);
   // namespace pt = boost::posix_time;
   // pt::ptime now = pt::second_clock::local_time();
    
   // std::string now_str = to_simple_string(now);
    
    std::cout << std::asctime(std::localtime(&t)) << "\tDEBUG\t" << m << std::endl;
    
}

void error(std::string m) {
    // formatting    : time [idx] message \n
    // destinations  : console, file "out.txt" and debug window
    
    std::time_t t = std::time(nullptr);
    
    //namespace pt = boost::posix_time;
    //pt::ptime now = pt::second_clock::local_time();
    
    //std::string now_str = to_simple_string(now);
    
    std::cout << std::asctime(std::localtime(&t)) << "\tERROR\t" << m << std::endl;
    
}

Test::Test(){
    t = 10;
    throw 1;
}

Test::~Test(){
    
}

