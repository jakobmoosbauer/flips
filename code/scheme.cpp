#include "scheme.hpp"

Scheme::Scheme(string filename, int n, int m, int l) : lastFlip(0,0,0){
  this->n = n;
  this->m = m;
  this->l = l;
  
  int maxrank = n*m*l;
  rows = new mat*[maxrank];
  rank = 0;

  string line;
  string a,b,c;
  ifstream input(filename);
  while(getline(input, line)){
    int d1 = line.find('*');
    int d2 = line.find('*',d1+1);
    a = line.substr(0,d1);
    b = line.substr(d1+1, d2-d1-1);
    c = line.substr(d2+1);
    rows[rank] = new mat[3];
    parseMatrix(a,'a',m,rows[rank]);
    parseMatrix(b,'b',l,rows[rank]+1);
    parseMatrix(c,'c',n,rows[rank]+2);
    ++rank;
  }
  input.close();
}

Scheme::~Scheme(){
  for(int i = 0; i<rank; ++i){
    delete rows[i];
  }
  delete rows;
}

Scheme::Scheme(const Scheme &scheme) : lastFlip(0,0,0){
  n = scheme.n;
  m = scheme.m;
  l = scheme.l;

  rank = scheme.rank;
  rows = new mat*[rank];
  
  for(int i = 0; i<rank; ++i){
    rows[i]=new mat[3];
    rows[i][0] = scheme.rows[i][0];
    rows[i][1] = scheme.rows[i][1];
    rows[i][2] = scheme.rows[i][2];
  }
}

void parseMatrix(string s, char x, int m, mat* matrix){
  *matrix = 0;
  int pos = -1;
  while((pos = s.find(x,pos+1))!=string::npos){
    set(*matrix,m,s[pos+1]-'0',s[pos+2]-'0');
  }
}

void Scheme::write(string filename){
  ofstream output(filename);
  for(auto r=0; r<rank; ++r){
    writeMatrix(output,'a',n,m,rows[r][0]);
    output << '*';
    writeMatrix(output,'b',m,l,rows[r][1]);
    output << '*';
    writeMatrix(output,'c',l,n,rows[r][2]);
    output << endl;
  }
  output.close();
}

void Scheme::writetoconsole(){
  ostringstream output;
  for(auto r=0; r<rank; ++r){
    writeMatrix(output,'a',n,m,rows[r][0]);
    output << '*';
    writeMatrix(output,'b',m,l,rows[r][1]);
    output << '*';
    writeMatrix(output,'c',l,n,rows[r][2]);
    output << endl;
  }
  output << endl;
  cout << output.str();
}

void writeMatrix(ostream &output, char x, int n, int m, mat matrix){
  output << '(';
  for(int i = 1; i <= n; ++i){
    for(int j = 1; j <= m; ++j){
      if(get(matrix, m, i, j)){
	output << x << i << j << '+';
      }
    }
  }
  output.seekp(-1,ios::cur);
  output << ')';
}

Flip::Flip(int col, int row1, int row2){
  this->col = col;
  this->row1 = row1;
  this->row2 = row2;
}

vector<Flip> Scheme::findFlips(){
  vector<Flip> out;
  for(auto k = 0; k < 3; ++k){
    for(auto i = 0; i < rank; ++i){
      for(auto j = i+1; j < rank; ++j){
	if(rows[i][k] == rows[j][k]){
	  #ifdef DEBUG
	  cout << "found flip, col:" << k << ", row1:" << i << ", row2:" << j << endl;
	  #endif
	  out.push_back(Flip(k,i,j));
	}
      }
    }
  }
  return out;
}

void Scheme::flip(Flip f){
  flip(f.col,f.row1,f.row2);
}

void Scheme::flip(int col, int row1, int row2){
#ifdef DEBUG
  cout << "flipping " << row1 << " and " << row2 << " over column " << col << endl;
#endif
  assume(0<=col && col<3);
  rows[row1][plus2mod3[col]] ^= rows[row2][plus2mod3[col]];
  rows[row2][plus1mod3[col]] ^= rows[row1][plus1mod3[col]];
}

void gauss(mat* matrix, int d){
  for(int j=0; j<d; ++j){
    mat pivot = matrix[j]&(-matrix[j]);
    for(int k=j+1; k<d+1; ++k){
      if(pivot&(matrix[k])){
	matrix[k] ^= matrix[j];
      }
    }
  }
}

void gauss2(mat* matrix, int d){
  for(int j=0; j<d; ++j){
    mat pivot = matrix[j]&(-matrix[j]);
    for(int k=j+1; k<d+2; ++k){
      if(pivot&(matrix[k])){
	matrix[k] ^= matrix[j];
      }
    }
  }
}

void gauss(mat* matrix, int* comb, int d){
  for(int j=0; j<d; ++j){
    mat pivot = matrix[j]&(-matrix[j]);
    for(int k=j+1; k<d+1; ++k){
      if(pivot&(matrix[k])){
	matrix[k] ^= matrix[j];
	comb[k] ^= comb[j];
      }
    }
  }
}

void fullgauss(mat* matrix, int* comb, int d){
  for(int j=0; j<d; ++j){
    mat pivot = matrix[j]&(-matrix[j]);
    for(int k=j+1; k<d; ++k){
      if(pivot&(matrix[k])){
	matrix[k] ^= matrix[j];
	comb[k] ^= comb[j];
      }
    }
  }
  for(int j=0; j<d; ++j){
    mat pivot = matrix[d-1-j]&(-matrix[d-1-j]);
    for(int k=j+1; k<d; ++k){
      if(pivot&(matrix[d-1-k])){
	matrix[d-1-k] ^= matrix[d-1-j];
	comb[d-1-k] ^= comb[d-1-j];
      }
    }
  }
}


string Scheme::newfilename(){
  mat s = 0;
  for(auto i = 0; i < rank; ++i){
    s+=rows[i][0]+rows[i][1]+rows[i][2];
    s<<=1;
    s%=9223372036854775807;
  }
  stringstream stream;
  stream << "j" << setfill('0') << setw(15) << hex << s << ".exp";
  string name = stream.str();
  return name;
}

void writelog(string oldfilename, string newfilename, int steps, int oldrank, int rank){
  ofstream log;
  log.open(newfilename.substr(0,3) + ".log", ios_base::app);
  log << time(NULL) << "," << oldfilename << "," << oldrank << "," << newfilename << "," << steps << "," << rank << endl;
  log.close();
}

//Test whether a scheme is a correct matrix multiplication scheme
bool Scheme::check(){
  mat** t = new mat*[64];
  for(int i = 0; i < n*m; ++i){
    t[i] = new mat[64];
    for(int j = 0; j < m*l; ++j){
      t[i][j] = 0;
    }
  }
  
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j){
      for(int k = 0; k < l; ++k){
	t[m*i+j][l*j+k] = (mat)1 << (n*k+i);
      }
    }
  }
  
  for(int s = 0; s < rank; ++s){
    for(int i = 0; i < n*m; ++i){
      for(int j = 0; j < m*l; ++j){
	if(rows[s][0]&(mat)1<<i && rows[s][1]&(mat)1<<j){
	  t[i][j] ^= rows[s][2];
	}
      }
    }
  }
  
  for(int i = 0; i < n*m; ++i){
    for(int j = 0; j < m*l; ++j){
      if(t[i][j] != 0){
	return false;
      }
    }
    delete[] t[i];
  }
  delete[] t;
  return true;
}



// reduce the scheme if posssible
bool Scheme::reduce(){
  unordered_multimap<mat,int>* xpos = new unordered_multimap<mat, int> [3];
  
  for(auto k = 0; k < rank; ++k){
    xpos[0].emplace(rows[k][0],k);
    xpos[1].emplace(rows[k][1],k);
    xpos[2].emplace(rows[k][2],k);
  }

  for(int a = 0; a<3; ++a){
    int b = plus1mod3[a];
    int c = plus1mod3[b];

    auto i = xpos[a].begin();
    while(i!= xpos[a].end()){
      auto range = xpos[a].equal_range(rows[i->second][a]);
      int d = distance(range.first, range.second);
      mat* bpool = new mat[d];
      mat* cpool = new mat[d];
      int* bcomb = new int[d];
      int* ccomb = new int[d];
      auto it = range.first;
      for(int j = 0; j<d; ++j, ++it){
	bpool[j] = rows[it->second][b];
	cpool[j] = rows[it->second][c];
	bcomb[j] = ccomb[j] = 1<<j;
      }
      fullgauss(bpool, bcomb, d);
      fullgauss(cpool, ccomb, d);
      for(int j = 0; j<d; ++j){
	if(bpool[j] == 0){
	  it = range.first;
	  for(int k = 0; k<j; ++k, ++it);
	  auto rj = rows[it->second];
	  it = range.first;
	  for(int k = 0; k<d; ++k, ++it){
	    if(bcomb[j]&(1<<k) && k != j){
	      rows[it->second][c] ^= rj[c];
	    }
	  }
	  for(it = range.first; j>0 ; --j, ++it);
	  rank--;
	  delete[] rows[it->second];
	  rows[it->second] = rows[rank];
	  return true;
	}
      }
      for(int j = 0; j<d; ++j){
	if(cpool[j] == 0){
	  it = range.first;
	  for(int k = 0; k<j; ++k, ++it);
	  auto rj = rows[it->second];
	  it = range.first;
	  for(int k = 0; k<d; ++k, ++it){
	    if(ccomb[j]&(1<<k) && k != j){
	      rows[it->second][b] ^= rj[b];
	    }
	  }
	  for(it = range.first; j>0 ; --j, ++it);
	  rank--;
	  delete[] rows[it->second];
	  rows[it->second] = rows[rank];
	  return true;
	}
      }
      i = range.second;
    }
  }
  return false;
}

void Scheme::randompathwithoutreduction(int steps){
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  srand((long)ts.tv_nsec);
  for(int i=0;i<steps;++i){
    vector<Flip> flips = findFlips();
    flip(flips[rand()%flips.size()]);
  }
  reduce();
}

void mylog(string s){
  ofstream errorlog;
  errorlog.open("errorlog.txt",ios_base::app);
  errorlog << s << endl;
  errorlog.close();
}

bool operator==(const Flip& f1, const Flip& f2){
  return(f1.row1 == f2.row1 && f1.row2 == f2.row2 && f1.col == f2.col);
}

