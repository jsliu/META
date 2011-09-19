#include "parameter.h"

Parameter::Parameter(int num_para, char* para[])
  :snptest(false), chr_included(false), sample_size(false), rsid(false), interval(false), top(false), gc(false), debug(false), zflag(false)
{

    /* read the command line values */
    map<string, int> para_exist;
    para_exist["--output"] = 0;
    para_exist["--cohort"] = 0;
    //para_exist["--sample-size"] = 0;
    para_exist["--method"] = 0;

    threshold = 0.5;    // default value

    for(int i = 1; i < num_para;) {
    	if (strcmp(para[i], "--output") == 0) {
    		para_exist[para[i]] = 1;
    		output = para[i+1];
		i += 2;
    	}
    	else if (strcmp(para[i], "--snptest") == 0) { snptest = true; i += 2; }
    	else if (strcmp(para[i],"--threshold") == 0) { threshold = string2double(para[i+1]); i += 2;}
    	else if (strcmp(para[i],"--method") == 0) {
    		para_exist[para[i]]= 1;
    		methods = string2int(para[i+1]);
    		if (methods == 3) {
    			for (int j = 1; j < num_para; j++) {
    				if (strcmp(para[j],"--sample-size") == 0) {
    					int k = 1;
    					while ((j + k) < num_para && *para[j + k] != '-') {
    						sample_sizes.push_back(string2long(para[j + k]));
    						k++;
    					}
    					sample_size = true;
					i += k;
    				}
    			}
    		}

    		if (methods > 3) {
    			string error = "The value of method must be in between 1 and 3.\n";
    			throw(BadFile(error));
    		}
		i += 2;
    	}
    	else if (strcmp(para[i],"--cohort") == 0) {
    		int j = 1;
    		while((i + j) < num_para && *para[i + j] != '-') {
    			cohorts.push_back(para[i + j]);
    			j++;
    		}
    		para_exist[para[i]] = 1;
		i += j;
    	}
    	else if (strcmp(para[i],"--rsid") == 0) {
    		int j = 1;
    		while((i + j) < num_para && *para[i + j] != '-') {
    			rsids.push_back(para[i + j]);
    			j++;
    		}
    		rsid = true;
		i += j;
    	}
    	else if (strcmp(para[i],"--interval") == 0) {
    		int j = 1;
    		while((i + j) < num_para && *para[i + j] != '-') {
    			region.push_back(string2long(para[i + j]));
    			j++;
    		}
#ifdef DEBUG
    		assert(region[0] <= region[1]);
#endif
    		if (region.size() != 2) {
    			string error = "Only two bouds are in need\n";
    			throw(BadFile(error));

    		}
    		interval = true;
		i += j;
    	}
    	else if (strcmp(para[i],"--top-snp") == 0) {
    		top_snp = string2long(para[i+1]);
    		top = true;
		i += 2;
    	}
    	else if (strcmp(para[i],"--lambda") == 0) {
    		int j = 1;
    		while ((i + j) < num_para && *para[i + j] != '-') {
    			genomic_controls.push_back(string2double(para[i + j]));
    			j++;
    		}
    		gc = true;
		i += j;
    	}
	else {
		string option(para[i]);
		string error = "ERROR! Parameter [ " + option + " ] is unrecognised.";
		throw BadFile(error);
	}
    }

    for(map<string, int>::iterator it_para = para_exist.begin(); it_para != para_exist.end(); it_para++) {

      if (!(it_para -> second)) {

        string error = "ERROR! Parameter [ " + it_para->first + " ] is in need.";
        throw BadFile(error);

      }

    }

    print();

}

void Parameter::print()
{

  cout << endl;
  cout << "         META v1.3.2 (beta)" << endl;
  cout << "===================================" << endl << endl;

  if (snptest) {
    cout << "==> SNPTEST output used" << endl;
  } else {
    cout << "==> Standard output used" << endl;
  }

  cout << "==> Method being used: " ;
  switch(methods) {
      case 1:
          cout << "inverse-variance method (fixed effects model)" << endl;
          break;
      case 2:
          cout << "inverse-variance method (random effects model)" << endl;
          break;
      case 3:
          cout << "combining z-statistics method" << endl;
          break;
  }

  if (rsid) {
    cout << "==> Selected SNPs: ";
    for(int i = 0; i < rsids.size(); i++)
      cout << "        " << rsids[i] << endl;
    cout << endl;
  }

  if (interval) {
    cout << "==> Selected Region: ";
    cout << " [ " << region[0] << " , " << region[1] << " ]  " << endl << endl;
  }
  cout << "==> OUTPUT: " << output << endl << endl;
}
