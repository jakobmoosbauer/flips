/***********************************************************************
pairSet.cpp

Copyright 2023 Jakob Moosbauer
 **********************************************************************/

#include "pairSet.hpp"

void PairSet::insert(uint64_t first, uint64_t second){
  uint64_t elem = (first << 32) | second;
  pairs.push_back(elem);
}

void PairSet::remove(uint64_t first, uint64_t second){
  uint64_t pair = (first << 32) | second;
  uint64_t pair2 = (second << 32) | first;
  for(size_t i = 0; i < pairs.size(); ++i){
    if(pairs[i] == pair || pairs[i] == pair2){
      pairs[i] = pairs.back();
      pairs.pop_back();
      return;
    }
  }
}

void PairSet::remove(uint64_t elem){
  for(size_t i = 0; i < pairs.size(); ++i){
    if(first(i) == elem || second(i) == elem) {
      pairs[i] = pairs.back();
      pairs.pop_back();
      --i;
    }
  }
}

size_t PairSet::size(){
  return pairs.size();
}

bool PairSet::contains(uint64_t first, uint64_t second){
  uint64_t pair = (first << 32) + second;
  uint64_t pair2 = (second << 32) + first;
  for(auto i :pairs){
    if(i == pair || i == pair2){
      return true;
    }
  }
  return false;
}

bool PairSet::contains(uint64_t pair){
  uint64_t pair2 = (pair & 0xFFFFFFFF) | (pair >> 32);
  for(auto i : pairs){
    if(i == pair || i == pair2){
      return true;
    }
  }
  return false;
}

void PairSet::clear(){
  pairs.clear();
}

uint64_t PairSet::first(size_t i){
  return pairs[i] >> 32;
}

uint64_t PairSet::second(size_t i){
  return pairs[i] & 0xFFFFFFFF;
}
