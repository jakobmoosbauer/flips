/***********************************************************************
mm.hpp

Copyright 2023 Jakob Moosbauer
 **********************************************************************/

#ifndef mm_hpp___
#define mm_hpp___

#include "tensor.hpp"

#define getm(matrix,m,i,j) matrix&(((factor)1)<<(m*(i-1)+j-1))
#define setm(matrix,m,i,j) matrix = matrix|(((factor)1)<<(m*(i-1)+j-1))
#define unset(matrix,m,i,j) matrix = matrix&~(((factor)1)<<(m*(i-1)+j-1))

extern int correctness_check;

class MM: public Tensor{
public:
  int n;
  int m;
  int l;

  MM(string filename, int n);
  MM(string filename, int n, int m, int l);
  MM(const MM &t);

  virtual MM* clone() const;

  virtual void write(string filename);
  virtual void writetoconsole();

  virtual bool iscorrect();
};

void parseMatrix(string s, char x, int m, factor* f);
void writeMatrix(ostream &output, char x, int n, int m, factor f);

#endif
