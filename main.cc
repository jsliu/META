#include "meta_analysis.h"

int main(int argc, char* argv[])
{
  
  int parameter_number = argc; 

  if(parameter_number == 1) {

      cout << "META v1.3.3 - Meta Analysis\n";
      cout <<  "(c) 2009-2011 Jason Liu, University of Oxford\n\n";
    
      cout << "Options:\n";
      cout << "        --cohort FILEs              A vector of formatted files\n";
      cout << "        --method NUMBER             Meta analysis methods:\n";
      cout << "                                        = 1, inverse-variance method (based on fixed effects model)\n";
      cout << "                                        = 2, inverse-variance method (based on random effects model)\n";
      cout << "                                        = 3, z-statistics combining method\n";
      cout << "        --ouput FILE                Output file\n";
      cout << "        [--snptest]                 Use output of SNPTEST as input files\n";
      cout << "        [--sample-size] NUMBERs     A vector of sample sizes of each cohort\n";
      cout << "        [--threshold] NUMBER        Threshold of imputaton quality score (0 - 1), default value = 0.5\n";
      cout << "        [--rsid] RSIDs              RSIDs of SNPs of interest\n";
      cout << "        [--interval] NUMBERs        Lower and upper bounds of the region of interest\n";
      cout << "        [--top-snp] NUMBER          Number of most significant SNPs in the output file\n\n";
      cout << "Contact : jsliu@stats.ox.ac.uk\n\n";   

  }else{
    try {
      Parameter p(argc, argv);    
      META ma(p);
      ma.meta_run();
    }catch(exception& e){
      cerr << "Exiting: " << e.what() << endl;
    }

  }

  return 0;

}
