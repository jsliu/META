#ifndef _PARAMETER_H_
#define _PARAMETER_H_ 

#include <cassert>
#include <cstring>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <string>

#include "tools/utility.h"
#include "error.h"

using namespace std;


using namespace std;

typedef unsigned long Positive;

class Parameter
{
public:
  Parameter(int, char**);
  ~Parameter() {};
  void print();

  string output;
  int methods;
  double threshold;
  vector<string> cohorts;
  vector<long> sample_sizes;
  vector<string> rsids;
  vector<Positive> region;
  long top_snp;
  vector<double> genomic_controls;

  bool snptest;
  bool chr_included;
  bool sample_size;
  bool rsid;
  bool interval;
  bool top;
  bool gc;
  bool debug;
  bool zflag;
};

#endif
