/***********************************************************************
pairSet.hpp

Copyright 2023 Jakob Moosbauer
 **********************************************************************/

#ifndef pairset_hpp__
#define pairset_hpp__

#include<vector>
#include<cstdint>

using namespace std;

class PairSet{
  public:
  vector<uint64_t> pairs;

  void insert(uint64_t first, uint64_t second);
  void remove(uint64_t first, uint64_t second);
  void remove(uint64_t elem);
  size_t size();
  void clear();
  uint64_t first(size_t i);
  uint64_t second(size_t i);

  bool contains(uint64_t first, uint64_t second);
  bool contains(uint64_t pair);
};
  
#endif
