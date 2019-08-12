//
//  Float64HashTable.cpp
//  TEToolkit_c++
//
//  Created by Ying Jin on 2/17/17.
//  Copyright (c) 2017 Ying Jin. All rights reserved.
//

#include "Float64HashTable.h"

#include <iostream>
#include <vector>
#include <string>
#include "stdlib.h"


//A hashtable taking 64bit float as key and 64bit float as value.
//Float64HashTable::Float64HashTable()
//{
//    this->table = kh_init(64);

//}

Float64HashTable::Float64HashTable(int size_hint)
{
    this->table = kh_init(64);
    
        if (size_hint)
        {
            kh_resize(64,this->table, size_hint);
        }
    
}

int Float64HashTable::size()
{
    return (this->table)->n_buckets;
}
Float64HashTable::~Float64HashTable() //__dealloc__(self):
{
    kh_destroy(64,this->table);
}

//Check if a given key is valid in hashtable.
bool Float64HashTable::has_key(khfloat64_t key)
{

    khiter_t k;
    k = kh_get(64,this->table, key);
    return k != (this->table)->n_buckets;
}

//Given a key, return a value.
khfloat64_t Float64HashTable::get_item(khfloat64_t key)
{
    khiter_t k;
    k = kh_get(64,this->table, key);
    if (k != table->n_buckets){
        return table->vals[k];
    }
    else{
        //raise KeyError(key);
        std::cout << "get_item key error " << key << std::endl;
        exit(EXIT_FAILURE);
    }
}


//Put a key-value pair to hashtable.
void Float64HashTable::set_item(khfloat64_t key, khfloat64_t val)
{

    khiter_t k;
    int ret = 0;

    k = kh_put(64,table, key, &ret);
    
    table->keys[k] = key;
    
    if (kh_exist(table, k))
    {
        table->vals[k] = val;
    }
    else{
        std::cout << "set_item key error " << key << std::endl;
        exit(EXIT_FAILURE);
        //raise KeyError(key);
    }
}

//increase value
void Float64HashTable::inc_item(khfloat64_t key, khfloat64_t val)
{
    
    khiter_t k;
    int ret = 0;
    
    k = kh_put(64,table, key, &ret);
    
    table->keys[k] = key;
    
    if (kh_exist(table, k))
    {
        table->vals[k] += val;
    }
    else{
        std::cout << "inc_item key error " << key << std::endl;
        exit(EXIT_FAILURE);
        //raise KeyError(key);
    }
}

//Take an array of keys in float64, and an array of values in
//float64 (of the same length), create key-value pairs in hashtable.
void Float64HashTable::map(std::vector<khfloat64_t> keys, std::vector<khfloat64_t> values)
{
    int ret = 0;
    khfloat64_t key;
    khiter_t k;

    for(size_t i =0; i < values.size() ; i++){
        key = keys[i];
        k = kh_put(64,table, key, &ret);
        table->keys[k] = key;
        table->vals[k] =  values[i];
    }
}


//Take a list of keys, return their values.
std::vector<khfloat64_t> Float64HashTable::lookup(std::vector<khfloat64_t> values)
{
    
    khfloat64_t key;
    khiter_t k;
    std::vector<khfloat64_t> locs ;
    
    for (size_t j = 0; j < values.size(); j++) {
        locs.push_back(0.0);
    }

    for(size_t i = 0 ; i < values.size(); i++)
    {
        key = values[i];
        k = kh_get(64,table, key);
        
        if( k != table->n_buckets){
            locs[i] = table->vals[k];
        }
        else{
            locs[i] = -1;
        }
    }

    return locs;
}
