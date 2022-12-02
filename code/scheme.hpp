/****************
    Copyright 2022 Jakob Moosbauer


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

******************/

#ifndef scheme_hpp___
#define scheme_hpp___

#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <iterator>
#include <unordered_set>
#include <random>
#include <cstdlib>
#include <stdlib.h>
#include <sstream>
#include <iomanip>

#define get(matrix,m,i,j) matrix&(((mat)1)<<(m*(i-1)+j-1))
#define set(matrix,m,i,j) matrix = matrix|(((mat)1)<<(m*(i-1)+j-1))
#define unset(matrix,m,i,j) matrix = matrix&~(((mat)1)<<(m*(i-1)+j-1))

#ifdef EXPECT
#define ife(expr) if(__builtin_expect(expr,false))
#define assume(expr) __builtin_assume(expr)
#else
#define ife(expr) if(expr)
#define assume(expr) do {} while(0)
#endif

#ifdef DEBUG2
#define DEBUG
#endif

#ifdef DEBUG3
#define DEBUG
#endif

#ifdef DEBUG
#define debug(msg) cout << msg << endl
#else
#define debug(msg) do {} while(0)
#endif

using namespace std;

typedef unsigned long long mat;

static const int plus1mod3[] = {1,2,0};
static const int plus2mod3[] = {2,0,1};

class Flip{
public:
  int col;
  int row1;
  int row2;
  
  Flip(int, int, int);
};

class Scheme{
public:
  int n;
  int m;
  int l;
  int rank;
  mat** rows;

  Flip lastFlip;

  Scheme(string filename, int n);
  Scheme(string filename, int n, int m, int l);
  Scheme(const Scheme &scheme);

  ~Scheme();

  void write(string filename);
  void writetoconsole();

  string newfilename();

  vector<Flip> findFlips();

  void flip(int col, int row1, int row2);
  void flip(Flip f);

  bool check(mat p, mat q, mat t, int b, int c, int a, int i, unordered_multimap<mat,int> cpos);
  bool reduce3rows(mat*, mat*, mat*, int, int , int);
  bool reduce4rows(mat*, mat*, mat*, mat*, mat, int, int, int, int);
  void reduce(mat*, mat*, unordered_multimap<mat,int>::iterator, mat, int, int, int, int, int j);

  bool reduce();
  
  int randompath(int steps, long seed);
  void randompathwithoutreduction(int steps);

  bool check();
};

void parseMatrix(string s, char x, int m, mat* matrix);
void writeMatrix(ostream &output, char x, int n, int m, mat matrix);
 
void writelog(string, string, int, int, int);

void gauss(mat*, int);
void gauss2(mat*, int);
void gauss(mat*, int*, int);

void mylog(string);

bool operator==(const Flip& f1, const Flip& f2);

#endif
