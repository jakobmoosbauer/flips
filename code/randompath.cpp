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

#include "scheme.hpp"

int Scheme::randompath(int steps, long seed){
  srand(seed);
  unordered_multimap<mat,int>* xpos = new unordered_multimap<mat, int> [3];
  
  for(auto k = 0; k < rank; ++k){
    xpos[0].emplace(rows[k][0],k);
    xpos[1].emplace(rows[k][1],k);
    xpos[2].emplace(rows[k][2],k);
  }

  //check for instant reduction
  if(reduce()){
    return 0;
  }

  vector<Flip> flips;
  for(auto step = 1; step <= steps; ++step){
    flips.clear();
    for(auto a = 0; a < 3; ++a){
      int b = plus1mod3[a];
      int c = plus1mod3[b];
      
      auto i = xpos[a].begin();
      while(i != xpos[a].end()){
	mat* row = rows[i->second];
	auto last = xpos[a].equal_range(row[a]).second;

	// compute all possible flips
	for(auto j = i; j != last; ++j){
	  auto k = j;
	  ++k;
	  mat* rj = rows[j->second];
	  for(; k != last; ++k){
	    mat* rk = rows[k->second];             // rj[a] == rk[a] is ensured

	    //check whether a flip leads to a reduction
	    
	    //check if rows matching rj[b]^rk[b] can be linearly combined
	    auto range = xpos[b].equal_range(rj[b]^rk[b]);
	    int d = distance(range.first,range.second);
	    assume(d>=0);
	    if(d!=0){
	      mat* ri1 = rows[range.first->second];
	      if(d==1){
		ife(rj[c] == ri1[c]){
		  flip(a, k->second, j->second);
		  if(!reduce()){
		    mylog("49");
		    return -5;
		  }
		  return step;
		}
		ife(rk[c] == ri1[c]){
		  flip(a, j->second, k->second);
		  if(!reduce()){
		    mylog("57");
		    return -5;
		  }
		  return step;
		}
	      }
	      else if(d==2){
		mat* ri2 = rows[(++range.first)->second];
		ife(rj[a] == (ri1[a]^ri2[a])){
		  flip(a,j->second,k->second);
		  if(!reduce()){
		    mylog("68");
		    return -5;
		  }
		  return step;
		}
		ife(rk[c] == (ri1[c]^ri2[c])){
		  flip(a,j->second,k->second);
		  if(!reduce()){
		    mylog("76");
		    return -5;
		  }
		  return step;
		}
		ife(rj[c] == (ri1[c]^ri2[c])){
		  flip(a,k->second,j->second);
		  if(!reduce()){
		    mylog("84");
		    return -5;
		  }
		  return step;
		}
		ife(rj[c] == ri1[c] || rj[c] == ri2[c]){
		  flip(a,k->second,j->second);
		  if(!reduce()){
		    mylog("92");
		    return -5;
		  }
		  return step;
		}
		ife(rk[c] == ri1[c] || rk[c] == ri2[c]){
		  flip(a,j->second,k->second);
		  if(!reduce()){
		    mylog("100");
		    return -5;
		  }
		  return step;
		}
	      }
	      else{
		mat* apool = new mat[d+1];
		mat* cpool = new mat[d+2];
		auto it = range.first;
		for(int l = 0; l<d; ++it, ++l){
		  apool[l] = rows[it->second][a];
		  cpool[l] = rows[it->second][c];
		}
		cpool[d] = rj[c];
		cpool[d+1] = rk[c];
		apool[d] = rj[a];
		gauss2(cpool, d);
		gauss(apool,d);
		ife(cpool[d]==0){
		  flip(a,k->second, j->second);
		  if(!reduce()){
		    mylog("122");
		    return -5;
		  }
		  return step;
		}
		ife(cpool[d+1]==0){
		  flip(a,j->second, k->second);
		  if(!reduce()){
		    mylog("130");
		    return -5;
		  }
		  return step;
		}
		ife(apool[d]==0){
		  flip(a,k->second, j->second);
		  if(!reduce()){
		    mylog("138");
		    return -5;
		  }
		  return step;
		}
		delete[] apool;
		delete[] cpool;
	      }
	    }

	    //check if rows matching rj[c]^rk[c] can be linearly combined
	    range = xpos[c].equal_range(rj[c]^rk[c]);
	    d = distance(range.first,range.second);
	    assume(d>=0);
	    if(d > 1){
	      if(d==2){
		mat* ri1 = rows[range.first->second];
		mat* ri2 = rows[(++range.first)->second];
		ife(rj[b] == (ri1[b]^ri2[b])){
		  flip(a,j->second,k->second);
		  if(!reduce()){
		    mylog("159");
		    return -5;
		  }
		  return step;
		}
		ife(rk[b] == (ri1[b]^ri2[b])){
		  flip(a,k->second,j->second);
		  if(!reduce()){
		    mylog("167");
		    return -5;
		  }
		  return step;
		}
	      }
	      else{
		mat* bpool = new mat[d+2];
		auto it = range.first;
		for(int l = 0; l<d; ++it, ++l){
		  bpool[l] = rows[it->second][b];
		}
		bpool[d] = rj[b];
		bpool[d+1] = rk[b];
		gauss2(bpool, d);
		ife(bpool[d]==0){
		  flip(a,j->second,k->second);
		  if(!reduce()){
		    mylog("185");
		    return -5;
		  }
		  return step;
		}
		ife(bpool[d+1]==0){
		  flip(a,k->second,j->second);
		  if(!reduce()){
		    mylog("193");
		    return -5;
		  }
		  return step;
		}
		delete[] bpool;
	      }
	    }

	    //check if rows matching rj[c] can be linearly combined
	    range = xpos[c].equal_range(rj[c]);
	    d = distance(range.first,range.second)-1;
	    assume(d>=0);
	    if(d>1){
	      if(d==2){
		mat* ri1;
		mat* ri2;
		if(range.first->second == j->second){
		  ri1 = rows[(++range.first)->second];
		  ri2 = rows[(++range.first)->second];
		}
		else{
		  ri1 = rows[range.first->second];
		  if((++range.first)->second == j->second){
		    ri2 = rows[(++range.first)->second];
		  }
		  else{
		    ri2 = rows[range.first->second];
		  }
		}
		ife((rj[b]^rk[b]) == (ri1[b]^ri2[b])){
		  flip(a,k->second,j->second);
		  if(!reduce()){
		    mylog("226");
		    return -5;
		  }
		  return step;
		}
	      }
	      else{
		mat* bpool = new mat[d+1];
		auto it = range.first;
		for(int l = 0; l<d; ++it, ++l){
		  if(it->second == j->second){
		    ++it;
		  }
		  bpool[l] = rows[it->second][b];
		}
		bpool[d] = rj[b]^rk[b];
		gauss(bpool,d);
		ife(bpool[d]==0){
		  flip(a,k->second,j->second);
		  if(!reduce()){
		    mylog("248");
		    return -5;
		  }
		  return step;
		}
		delete[] bpool;
	      }
	    }

	    //check if rows matching rk[c] can be linearly combined
	    range = xpos[c].equal_range(rk[c]);
	    d = distance(range.first,range.second)-1;
	    assume(d>=0);
	    if(d>1){
	      if(d==2){
		mat* ri1;
		mat* ri2;
		if(range.first->second == k->second){
		  ri1 = rows[(++range.first)->second];
		  ri2 = rows[(++range.first)->second];
		}
		else{
		  ri1 = rows[range.first->second];
		  if((++range.first)->second == k->second){
		    ri2 = rows[(++range.first)->second];
		  }
		  else{
		    ri2 = rows[range.first->second];
		  }
		}
		ife((rj[b]^rk[b]) == (ri1[b]^ri2[b])){
		  flip(a,j->second,k->second);
		  if(!reduce()){
		    mylog("279");
		    return -5;
		  }
		  return step;
		}
	      }
	      else{
		mat* bpool = new mat[d+1];
		auto it = range.first;
		for(int l = 0; l<d; ++it, ++l){
		  if(it->second == k->second){
		    ++it;
		  }
		  bpool[l] = rows[it->second][b];
		}
		bpool[d] = rj[b]^rk[b];
		gauss(bpool,d);
		ife(bpool[d]==0){
		  flip(a,j->second,k->second);
		  if(!reduce()){
		    mylog("299");
		    return -5;
		  }
		  return step;
		}
		delete[] bpool;
	      }
	    }

	    //check if rows matching rj[b] can be linearly combined
	    range = xpos[b].equal_range(rj[b]);
	    d = distance(range.first,range.second)-1;
	    assume(d>=0);
	    if(d>1){
	      if(d==2){
		mat* ri1;
		mat* ri2;
		if(range.first->second == j->second){
		  ri1 = rows[(++range.first)->second];
		  ri2 = rows[(++range.first)->second];
		}
		else{
		  ri1 = rows[range.first->second];
		  if((++range.first)->second == j->second){
		    ri2 = rows[(++range.first)->second];
		  }
		  else{
		    ri2 = rows[range.first->second];
		  }
		}
		ife((rj[c]^rk[c]) == (ri1[c]^ri2[c])){
		  flip(a,j->second,k->second);
		  if(!reduce()){
		    mylog("332");
		    return -5;
		  }
		  return step;
		}
	      }
	      else{
		mat* cpool = new mat[d+1];
		auto it = range.first;
		for(int l = 0; l<d; ++it, ++l){
		  if(it->second == j->second){
		    ++it;
		  }
		  cpool[l] = rows[it->second][c];
		}
		cpool[d] = rj[c]^rk[c];
		gauss(cpool,d);
		ife(cpool[d]==0){
		  flip(a,j->second,k->second);
		  if(!reduce()){
		    mylog("352");
		    return -5;
		  }
		  return step;
		}
		delete[] cpool;
	      }
	    }

	    //check if rows matching rk[b] can be linearly combined
	    range = xpos[b].equal_range(rk[b]);
	    d = distance(range.first,range.second)-1;
	    assume(d>=0);
	    if(d>1){
	      if(d==2){
		mat* ri1;
		mat* ri2;
		if(range.first->second == k->second){
		  ri1 = rows[(++range.first)->second];
		  ri2 = rows[(++range.first)->second];
		}
		else{
		  ri1 = rows[range.first->second];
		  if((++range.first)->second == k->second){
		    ri2 = rows[(++range.first)->second];
		  }
		  else{
		    ri2 = rows[range.first->second];
		  }
		}
		ife((rj[c]^rk[c]) == (ri1[c]^ri2[c])){
		  flip(a,k->second,j->second);
		  if(!reduce()){
		    mylog("385");
		    return -5;
		  }
		  return step;
		}
	      }
	      else{
		mat* cpool = new mat[d+1];
		auto it = range.first;
		for(int l = 0; l<d; ++it, ++l){
		  if(it->second == k->second){
		    ++it;
		  }
		  cpool[l] = rows[it->second][c];
		}
		cpool[d] = rj[c]^rk[c];
		gauss(cpool,d);
		ife(cpool[d]==0){
		  flip(a,k->second,j->second);
		  if(!reduce()){
		    mylog("405");
		    return -5;
		  }
		  return step;
		}
		delete[] cpool;
	      }
	    }
	    
	    flips.push_back(Flip(a,j->second, k->second));
	    flips.push_back(Flip(a,k->second, j->second));
	  }
	}
	i = last;
      }
    }
    if(flips.size() == 0){
      return -1;
    }

    #ifdef CONTROL
    for(Flip f:flips){
      flip(f);
      if(reduce())
	return -6;
      flip(f);
    }
    #endif

    Flip f = flips[rand()%flips.size()];
    while(f == lastFlip){
      f = flips[rand()%flips.size()];
    }
    lastFlip = f;
    int a = f.col;
    int b = plus1mod3[f.col];
    int c = plus2mod3[f.col];
    auto brange = xpos[b].equal_range(rows[f.row2][b]);
    auto i = brange.first;
    while(i->second != f.row2){
      ++i;
    }
    xpos[b].erase(i);
    auto crange = xpos[c].equal_range(rows[f.row1][c]);
    auto j = crange.first;
    while(j->second != f.row1){
      ++j;
    }
    xpos[c].erase(j);
    flip(f);
    xpos[b].emplace(rows[f.row2][b],f.row2);
    xpos[c].emplace(rows[f.row1][c],f.row1);
  }

  return -1;
}
