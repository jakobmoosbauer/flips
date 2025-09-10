/***********************************************************************
tensor.hpp

Copyright 2023 Jakob Moosbauer
 **********************************************************************/

#ifndef tensor_hpp___
#define tensor_hpp___

#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <random>
#include <cstdlib>
#include <stdlib.h>
#include <sstream>
#include "pairSet.hpp"
#include <iomanip>
#include <csignal>

#ifdef DEBUG
#define debug(msg) cerr << msg << endl
#else
#define debug(msg) do {} while(0)
#endif

using namespace std;

typedef unsigned long long factor;

static const int plus1mod3[] = {1,2,0};
static const int plus2mod3[] = {2,0,1};

extern int oldrank;
extern string filename;
extern int correctness_check;
//extern volatile sig_atomic_t termination_flag;

class Tensor{
public:
  int rank;
  int maxrank;
  factor* data;
  PairSet* flips;

  Tensor();
  Tensor(const Tensor &t);
  
  ~Tensor();

  virtual Tensor* clone() const;
  
  virtual void write(string filename);
  virtual void writetoconsole();
  
  virtual string newfilename();
  void writetofile(int);

  factor& get(int, int);
  void remove(int);
  void init();

  bool flip(int col, int row1, int row2, bool reduce_flag = true);
  void split(int col, int row1, int row2);
  bool reduce();
  void remove_zero_rows();

  bool randomflip(mt19937 &gen, uniform_int_distribution<> &coinflip, bool reduce_flag = true);
  
  void randompath(int steps, mt19937 &gen, int split_distance, bool split, bool restart);
  bool randomsplit(mt19937 &gen, uniform_int_distribution<> &coinflip, uniform_int_distribution<> &d3, int split_distance);

  virtual bool iscorrect();
};

void writelog(string, string, int, int);



#endif
