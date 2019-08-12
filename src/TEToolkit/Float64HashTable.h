//
//  Float64HashTable.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/17/17.
//  Copyright (c) 2017 Ying Jin. All rights reserved.
//

#ifndef __TEToolkit_c____Float64HashTable__
#define __TEToolkit_c____Float64HashTable__

#include <stdio.h>
#include <vector>
#include "khash.h"


KHASH_MAP_INIT_FLOAT64(64, khfloat64_t)

class Float64HashTable{
public:
    khash_t(64) *table;
    
    //Float64HashTable();
    Float64HashTable(int size_hint=1);
    ~Float64HashTable();
    
    std::vector<khfloat64_t> lookup(std::vector<khfloat64_t> values);
    void map(std::vector<khfloat64_t> keys, std::vector<khfloat64_t> values);
    void set_item(khfloat64_t key, khfloat64_t val);
    void inc_item(khfloat64_t key, khfloat64_t val);
    khfloat64_t get_item(khfloat64_t key);
    bool has_key(khfloat64_t key);
    int size();
    
};
#endif /* defined(__TEToolkit_c____Float64HashTable__) */
