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

# include "scheme.hpp"

int main(int argc, char* argv[]){
  debug("debugging enabled");

  int n,m,l;
  int pathlength;
  int tries;
  
  if(argc == 5){
    n=m=l=(*argv[2]-'0');
    pathlength = strtol(argv[3],NULL,10);
    tries = strtol(argv[4],NULL,10);
  }
  else if(argc == 7){
    n=*argv[2]-'0';
    m=*argv[3]-'0';
    l=*argv[4]-'0';
    pathlength = strtol(argv[5],NULL,10);
    tries = strtol(argv[6],NULL,10);
  }
  else{
    cerr << "Wrong number of arguments. Usage: " << endl;
    cerr << argv[0] << " <filename> <matrix size> <path length> <retries> or" << endl;
    cerr << argv[0] << " <filename> <dim 1> <dim 2> <dim 3> <path length> <retries>" << endl;
    return 1;
  }

  string filename = argv[1];

  Scheme s = Scheme(filename,n,m,l);

  int oldrank = s.rank;
  
  if(!s.check()){
    cerr << "Opened incorrect scheme: " << filename << endl;
    return 1;
  }
  
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  long seed = (long)ts.tv_nsec;

  
  int steps = -1;
  int stepsum = 0;

  if(tries == 0){
    stepsum = s.randompath(pathlength,seed);
  }
  else{
    do {
      for(int i = 0; i < tries; ++i){
	steps = s.randompath(pathlength, seed);
	if(steps >= 0){
	  stepsum += steps;
	  break;
	}
      }
    }
    while(steps >= 0);
  }

#ifdef DEBUG
  if(!s.check()){
    string newfilename = s.newfilename();
    s.write("errors/" + newfilename);
    cerr << "Found incorrect scheme. Incorrect scheme name written to errorlog.txt" << endl;
    mylog(filename + "," + "errors/" + newfilename + "," + to_string(steps) + ", seed:" + to_string(seed));
    return 1;
  }
#endif
  
  if(steps==-5){
    cerr << "Reduction not possible, scheme was not saved. Information written to errorlog.txt" << endl;
    return 1;
  }

  #ifdef CONTROL
  if(steps==-6){
    mylog("Found undetected reduction.");
  }
  #endif
  
  if(stepsum > 0){
    string newfilename = s.newfilename();
    s.write(newfilename);
    cout << newfilename << "," << s.rank << endl;
    writelog(filename, newfilename, stepsum, oldrank, s.rank);
  }
  return 0;
}
