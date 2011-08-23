// body file: meta_analysis.cc

#ifndef _META_ANALYSIS_
#define _META_ANALYSIS_


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <map>
#include <vector>
#include <algorithm>
#include <zlib.h>

#include "tools/cprob/cprob.h"
#include "tools/utility.h"
#include "error.h"

using namespace std;

extern "C" {

  double chdtrc(double v, double x);
  double ndtr(double x);
  double ndtri(double x);

}

typedef unsigned long int Positive;

class SNP
{

public:
  SNP() : chr(0), sample_size(0.), coded_af(0.) {};
  ~SNP() {};
  string rsid;
  char allele_A, allele_B;
  unsigned int chr;
  Positive pos;
  int imputed;
  long sample_size;
  double beta, se;
  double pvalue;
  double coded_af;

};


class CohortSNP : public SNP
{

public:
  CohortSNP() : genotype_aa(0.0), genotype_ab(0.0), genotype_bb(0.0) {};
  ~CohortSNP() {};
  double measure_info;
  double genotype_aa, genotype_ab, genotype_bb;
  bool snp_flipped;
  bool snp_reversed;

};


class MetaSNP : public SNP
{

public:
  MetaSNP() {};
  ~MetaSNP() {};
  double nu; // between study variance
  double Q;
  double p_het;
  double I2;

};


// Parameter for SNPTEST output
typedef map<Positive, CohortSNP> Cohort;
typedef map<Positive, MetaSNP> Meta;
typedef vector<Positive> Union;


class Parameter
{
public:
  Parameter(int, vector<string> &);
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
};

class META : public Parameter
{

public:
  META(int, vector<string> &);
  ~META() {};

  void read_data();
  gzFile read_cohort(string, bool quiet = true);
  void read_header(vector<string>&, const gzFile);
  void read_specific_snp(Cohort&, const CohortSNP&);
  void read_SNP_in_region(Cohort&, const CohortSNP&);
  int read_line(vector<string>&, const gzFile);
  int read_snptest_output(Cohort&, const gzFile);
  int read_standard_output(Cohort&, const gzFile);
  void select_SNP();
  void filter_data(Cohort&, const CohortSNP&);
  void printINFO(const int, const int, const int, const string);
  void print_missing_SNP();
  void write_data();
  void write_cohort(ofstream&);
  void write_SNP(ofstream&);
  void write_top_SNP(ofstream&);
  void meta_analysis();
  void meta_run();
  void create_union_list();
  void count_alleles();
  void check_alleles_across_cohorts(int&, int&, int&);
  void inverse_variance_method();
  void combine_z_score();
  void combine_beta_se_with_nu(const Positive, double&, double&);
  void combine_beta_se(const Positive, double&, double&);
  bool check_alleles_flipped(char, char, char, char);
  bool check_alleles_reversed(char, char, char, char);
  void heterogeneity(MetaSNP&, const double);

private:
  Meta meta;
  Meta::iterator it_meta;
  Union union_pos;
  Union::iterator it_union;

  MetaSNP meta_snp;
  //CohortSNP cohort_snp;

  vector<Cohort> cohort_list;
  vector<vector<double> > pvalue_list;
  vector<double> beta_i;
  vector<double> wei_i;
  vector<int> included; // to judge if the SNP of cohort is included at the position

  int number_cohort;
};


#endif
