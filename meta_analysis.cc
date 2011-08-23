#include "meta_analysis.h"

Parameter::Parameter(int num_para, vector<string>& para)
  :snptest(false), chr_included(false), sample_size(false), rsid(false), interval(false), top(false), gc(false)
{

    /* read the command line values */
    map<string, int> para_exist;
    para_exist["--output"] = 0;
    para_exist["--cohort"] = 0;
    //para_exist["--sample-size"] = 0;
    para_exist["--method"] = 0;

    threshold = 0.5;    // default value

    for(int i = 1; i < num_para; i++) {
    	if (para[i].find("--output") != string::npos) {
    		output = para[i + 1];
    		para_exist[para[i]] = 1;
    	}
    	if (para[i].find("--snptest") != string::npos) snptest = true;
    	if (para[i].find("--threshold") != string::npos) threshold = string2double(para[i+1]);
    	if (para[i].find("--method") != string::npos) {
    		methods = string2int(para[i + 1]);
    		if (methods == 3) {
    			for (int j = 1; j < num_para; j++) {
    				if (para[j].find("--sample-size") != string::npos) {
    					int k = 1;
    					while ((j + k) < num_para && para[j + k].find_first_of("-") == string::npos) {
    						sample_sizes.push_back(string2long(para[j + k]));
    						k++;
    					}
    					sample_size = true;
    				}
    			}
    		}

    		if (methods > 3) {
    			string error = "The value of method must be in between 1 and 3.\n";
    			throw(BadFile(error));
    		}
    		para_exist[para[i]]= 1;
    	}
    	if (para[i].find("--cohort") != string::npos) {
    		int j = 1;
    		while((i + j) < num_para && para[i + j].find_first_of("-") == string::npos) {
    			cohorts.push_back(para[i + j]);
    			j++;
    		}
    		para_exist[para[i]] = 1;
    	}
    	if (para[i].find("--rsid") != string::npos) {
    		int j = 1;
    		while((i + j) < num_para && para[i + j].find_first_of("-") == string::npos) {
    			rsids.push_back(para[i + j]);
    			j++;
    		}
    		rsid = true;
    	}
    	if (para[i].find("--interval") != string::npos) {
    		int j = 1;
    		while((i + j) < num_para && para[i + j].find_first_of("-") == string::npos) {
    			region.push_back(string2long(para[i + j]));
    			j++;
    		}
    		assert(region[0] <= region[1]);
    		if (region.size() != 2) {
    			string error = "Only two bouds are in need\n";
    			throw(BadFile(error));

    		}
    		interval = true;
    	}
    	if (para[i].find("--top-snp") != string::npos) {
    		top_snp = string2long(para[i + 1]);
    		top = true;
    	}
    	if (para[i].find("--lambda") != string::npos) {
    		int j = 1;
    		while ((i + j) < num_para && para[i + j].find_first_of("-") == string::npos) {
    			genomic_controls.push_back(string2double(para[i + j]));
    			j++;
    		}
    		gc = true;
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


META::META(int np, vector<string>& p)
  :Parameter(np, p)
{

  number_cohort = cohorts.size();
  #ifdef DEBUG 
	cout << "number_cohort = " << number_cohort << endl;
  #endif

}

/*******************************************************
 *  Open the output of SNPTEST and read it line by line
 *******************************************************/
void META::read_data()
{
	cout << "Reading the data ..." << endl << flush;
	cout << "==> " << number_cohort << " cohorts being used: " << endl << flush;

	/* to judge if the selected SNPs are in the data */
	if(rsid) select_SNP();

	/* read the data from all cohorts */
	cohort_list.reserve(number_cohort);
	for(vector<string>::iterator  it = cohorts.begin(); it != cohorts.end(); it++)
	{
		Cohort single_cohort;
		gzFile in = read_cohort(*it);

		/* read content */
		int num_row;
		if (snptest) {
			num_row = read_snptest_output(single_cohort, in);
		} else {
			num_row = read_standard_output(single_cohort, in);
		}

		int indx = static_cast<int> (it - cohorts.begin() + 1);
		int cohort_size = single_cohort.size();
		string cohort_name = *it;
		printINFO(num_row, indx, cohort_size, cohort_name);

		cohort_list.push_back(single_cohort);
	}

	if(rsid) print_missing_SNP();

}


gzFile META::read_cohort(string cohort, bool quiet)
{
	if (!quiet) cout << "Reading cohort: " + cohort << endl;

	gzFile gz_in = gzopen(cohort.c_str(), "rb");

	if (gz_in == NULL) {
		string err = "ERROR! Fail to open file: " + cohort;
		throw BadFile(err);
	}
	return gz_in;
}

int META::read_line(vector<string>& line, const gzFile in)
{
	int mem = 1048576;
	char* buf = new char[mem]; // 1Mb buffer
	char* str = gzgets(in, buf, mem);

	if (str == NULL) {
		delete []buf;
		return 0;
	}

	istringstream iss(str, istringstream::in);
	line.clear();
	while (!iss.eof()) {
		string strg;
		iss >> strg;
		line.push_back(strg);
	}
	delete []buf;
	return line.size() - 1;  // exclude the null terminator
}

int META::read_snptest_output(Cohort& single_cohort, const gzFile in)
{
	/* read header */
	vector<string> header;
	int header_column = read_line(header, in);

	int num_row = 0;
	vector<string> line;
	int body_column = read_line(line, in);
	while (body_column) {

		num_row++;
		if (body_column != header_column) {
			string err = "ERROR! The columns are not consistent with header columns at row " + int2string(num_row);
			throw BadFile(err);
		}

		CohortSNP cohort_snp;
		for(int i = 0; i < body_column; i++)
		{
			if(header[i] == "rsid") cohort_snp.rsid = line[i];
			if(header[i] == "pos") cohort_snp.pos = string2long(line[i]);
			if(header[i] == "allele_A") cohort_snp.allele_A = string2char(line[i]);
			if(header[i] == "allele_B") cohort_snp.allele_B = string2char(line[i]);
			if(header[i] == "all_AA") cohort_snp.genotype_aa = string2double(line[i]);
			if(header[i] == "all_AB") cohort_snp.genotype_ab = string2double(line[i]);
			if(header[i] == "all_BB") cohort_snp.genotype_bb = string2double(line[i]);
			
			if(header[i].find("_pvalue") != string::npos) cohort_snp.pvalue = string2double(line[i]);
			if(header[i].find("_info") != string::npos) cohort_snp.measure_info = string2double(line[i]);
			if(header[i].find("_beta_1") != string::npos) cohort_snp.beta = string2double(line[i]);
			if(header[i].find("_se_1") != string::npos) cohort_snp.se = string2double(line[i]);
		}
		//TODO: have to deal with the case when these values are not read
		filter_data(single_cohort, cohort_snp);
		body_column = read_line(line, in);
	} // while

	return num_row;
}

/* only read SNP with measure info > threshold */
int META::read_standard_output(Cohort& single_cohort, const gzFile in)
{
    /* read header, the 2, 3, 4, 5 and 6 column must be rsid, position, allele_A, allele_B and P_value */
	vector<string> header;
	int header_column = read_line(header, in);

	vector<string> line;
	int body_column = read_line(line, in);
	int num_row = 0;
	while (body_column) {

		num_row++;
		if (body_column != header_column) {
				string err = "ERROR! The columns are not consistent with header columns at " + int2string(num_row);
				throw BadFile(err);
		}

		CohortSNP cohort_snp;

    	for(int i = 0; i < body_column; i++)
    	{
    		if(header[i] == "chr") { cohort_snp.chr = string2int(line[i]); chr_included = true; }
    		if(header[i] == "rsid") cohort_snp.rsid = line[i];
    		if(header[i] == "pos") cohort_snp.pos = string2long(line[i]) + cohort_snp.chr * 1e9;
    		if(header[i] == "allele_A") cohort_snp.allele_A = string2char(line[i]);
    		if(header[i] == "allele_B") cohort_snp.allele_B = string2char(line[i]);
    		if(header[i] == "info") cohort_snp.measure_info = string2double(line[i]);
    		if(header[i] == "P_value") cohort_snp.pvalue = string2double(line[i]);
    		if(header[i] == "beta") cohort_snp.beta = string2double(line[i]);
    		if(header[i] == "se") cohort_snp.se = string2double(line[i]);
    	}
    	filter_data(single_cohort, cohort_snp);
    	body_column = read_line(line, in);
    } // while

	return num_row;
}


void META::select_SNP()
{
	    included.reserve(rsids.size());
	    for(int i = 0; i < rsids.size(); i++) included.push_back(1);

}


void META::filter_data(Cohort& single_cohort, const CohortSNP& cohort_snp)
{

	if(cohort_snp.measure_info >= threshold || cohort_snp.measure_info == -1)
	{

		if(cohort_snp.se > 0)
		{
			/* data for selected SNPs */
			if(rsid)
			{

				read_specific_snp(single_cohort, cohort_snp);

			}else if(interval) {

				read_SNP_in_region(single_cohort, cohort_snp);

			}else{

				single_cohort.insert(pair<Positive, CohortSNP>(cohort_snp.pos, cohort_snp));

			}

		}

    }

}

void META::read_specific_snp(Cohort& single_cohort, const CohortSNP& cohort_snp)
{

	for(int i = 0; i < rsids.size(); i++)
	{
		if(rsids[i].compare(cohort_snp.rsid) == 0)
		{
			single_cohort.insert(pair<Positive, CohortSNP>(cohort_snp.pos, cohort_snp));
			included[i]++;
		}
	}

}

void META::read_SNP_in_region(Cohort& single_cohort, const CohortSNP& cohort_snp)
{

	if(cohort_snp.pos >= region[0] && cohort_snp.pos <= region[1])
		single_cohort.insert(pair<Positive, CohortSNP>(cohort_snp.pos, cohort_snp));

}


void META::printINFO(const int num_row, const int indx, const int cohort_size, const string cohort_name)
{

	double percent = static_cast<double> (cohort_size) / static_cast<double> (num_row) * 100.0;

	string space = ": ";
	cout << "Cohort " << indx << space << cohort_name
		 << ", " << cohort_size << " out of " << num_row << " (";

	cout.precision(4);
	cout << percent;

	cout << "%) " << "SNPs are used (threshold >= " << threshold << ")" << endl;

}


void META::print_missing_SNP()
{

	for(int i = 0; i < rsids.size(); i++)
	{
		if(included[i] == 0)
		{

			string error = "OOPS! " + rsids[i] + " is not contained in all cohorts.";
	    	throw BadFile(error);

		}

	}

}
/********************************************
 *
 *   Print out the result of meta analysis
 *
 ********************************************/

void META::write_data()
{

	cout << "Writing the data ..." << endl << flush;

	ofstream out;
	write_cohort(out);

	if(chr_included) out << "chr" << " ";
	out << "rsid" << " ";
	out << "pos" << " ";
	out << "allele_A" << " ";
	out << "allele_B" << " ";
	out << "P_value" << " ";

	if(methods == 1 || methods == 2) {

		if(snptest) out << "coded_af" << " ";
		out << "beta" << " ";
		out << "se" << " ";
		out << "Q" << " ";
		out << "P_heterogeneity" << " ";
		out << "I2" << " ";

	}


	for(int j = 1; j <= number_cohort; j++) out << "P_cohort_" << j << " ";

	out << endl;

	/* write the output */
	#ifdef DEBUG
		cout << "meta size = " << meta.size() << endl;
	#endif

	if(top) {

		write_top_SNP(out);

	}else{

		write_SNP(out);

	}

	cout << endl << "DONE!" << endl << endl << flush;

}

void META::write_cohort(ofstream& out)
{
	out.open(output.c_str());

	if (out.fail()) {
	    string err = "ERROR! Fail to open file: " + output;
	    throw BadFile(err);
	}
}


void META::write_SNP(ofstream& out)
{
	int k = 0;

	for(it_meta = meta.begin(); it_meta != meta.end(); it_meta++) {

		if(chr_included) out << (*it_meta).second.chr << " ";
		out << (*it_meta).second.rsid << " ";
		out << (*it_meta).first << " ";
		out << (*it_meta).second.allele_A << " ";
		out << (*it_meta).second.allele_B << " ";
		out << (*it_meta).second.pvalue << " ";

		if(methods == 1 || methods == 2) {

			if(snptest) out << (*it_meta).second.coded_af << " ";
			out << (*it_meta).second.beta << " ";
			out << (*it_meta).second.se << " ";
			out << (*it_meta).second.Q << " ";
			out << (*it_meta).second.p_het << " ";
			out << (*it_meta).second.I2 << " ";

		}

		for(vector<double>::iterator it = pvalue_list[k].begin(); it != pvalue_list[k].end(); it++) out << *it << " ";

		out << endl;

		k++;

	}
}


void META::write_top_SNP(ofstream& out)
{
	map<double, MetaSNP, less<double> > sorted_snp;
	map<double, MetaSNP, less<double> >::iterator it_sorted;

	for(it_meta = meta.begin(); it_meta != meta.end(); it_meta++)
		sorted_snp[(*it_meta).second.pvalue] = (*it_meta).second;

	int k = 0;

	for(it_sorted = sorted_snp.begin(); it_sorted != sorted_snp.end(); it_sorted++) {

		if(k < top) {

			if(chr_included) out << (*it_sorted).second.chr << " ";

			out << (*it_sorted).second.rsid << " ";
			out << (*it_sorted).second.pos << " ";
			out << (*it_sorted).second.allele_A << " ";
			out << (*it_sorted).second.allele_B << " ";
			out << (*it_sorted).second.pvalue << " ";


			if(methods == 1 || methods == 2) {

				if(snptest) out << (*it_sorted).second.coded_af << " ";
				out << (*it_sorted).second.beta << " ";
				out << (*it_sorted).second.se << " ";
				out << (*it_sorted).second.Q << " ";
				out << (*it_sorted).second.p_het << " ";
				out << (*it_sorted).second.I2 << " ";

			}

			for(vector<double>::iterator it = pvalue_list[k].begin(); it != pvalue_list[k].end(); it++) out << *it << " ";

			out << endl;

		}

		k++;

	}


}


/************************************************
 *
 *   Meta analysis (by estimating effect size)
 *
 ************************************************/

void META::meta_analysis()
{

  create_union_list();

  if(methods == 1 | methods == 2) inverse_variance_method();
  if(methods == 3) combine_z_score();

}


void META::meta_run()
{

    read_data();
    meta_analysis();
    write_data();

}

/* create a union of snps in each cohort */
void META::create_union_list()
{

	const Positive bound = 1e+7;

	/* find the union of poss over all cohorts */
	Cohort::iterator it_cohort;
	union_pos.reserve(cohort_list[0].size());

	for(it_cohort = cohort_list[0].begin(); it_cohort != cohort_list[0].end(); it_cohort++)
	{

		#ifdef DEBUG
			cout << "cohort 1" << " ";
			cout << "rsid =  " << (*it_cohort).second.rsid << " ";
			cout << "pos = " << (*it_cohort).first << endl;
		#endif

		union_pos.push_back((*it_cohort).first);

	}


	if(number_cohort > 1) {

		vector<Positive> tmp_pos(bound);
		vector<Positive>::iterator vecit, tmp_it;

		for(int i = 1; i < number_cohort; i++) {

			vector<Positive> second_pos;
			second_pos.reserve(cohort_list[i].size());
			for(it_cohort = cohort_list[i].begin(); it_cohort != cohort_list[i].end(); it_cohort++) {

				#ifdef DEBUG
					cout << "cohort " << i + 1 << " ";
					cout << "rsid =  " << (*it_cohort).second.rsid << " ";
					cout << "pos = " << (*it_cohort).first << endl;
				#endif

				second_pos.push_back((*it_cohort).first);

			}

			vecit = set_union(union_pos.begin(), union_pos.end(), second_pos.begin(), second_pos.end(), tmp_pos.begin());

			union_pos.clear();
			union_pos.reserve(static_cast<int> (vecit - tmp_pos.begin()));
			for(tmp_it = tmp_pos.begin(); tmp_it != vecit; tmp_it++) union_pos.push_back(*tmp_it);

		}

	}

	cout << endl << "There are " << union_pos.size() << " SNPs in the union list." << endl;
	count_alleles();

}


/* check the alleles of snp for each cohort */
void META::count_alleles()
{
	int snp_flipped = 0;
	int snp_reversed = 0;
	int snp_removed = 0;

	for(it_union = union_pos.begin(); it_union < union_pos.end(); it_union++)
	{

		check_alleles_across_cohorts(snp_flipped, snp_reversed, snp_removed);

	}

	cout << snp_flipped << " SNPs are flipped!" << endl;
	cout << snp_reversed << " SNPs are reversed!" << endl;
	cout << snp_removed << " SNPs are removed!" << endl << endl;

}


void META::check_alleles_across_cohorts(int& snp_flipped, int& snp_reversed, int& snp_removed)
{
	char first_allele_A, first_allele_B, curr_allele_A, curr_allele_B;

	map<Positive, CohortSNP, less<Positive> >::iterator it_location;

	int first = 0;
	for(int i = 0; i < number_cohort; i++) {

		it_location = cohort_list[i].find(*it_union);
		if(it_location != cohort_list[i].end()) {

			if(first == 0) {

				first_allele_A = (*it_location).second.allele_A;
				first_allele_B = (*it_location).second.allele_B;

				(*it_location).second.snp_reversed = false;

			}else{

				curr_allele_A = (*it_location).second.allele_A;
				curr_allele_B = (*it_location).second.allele_B;

				(*it_location).second.snp_flipped = false;
				(*it_location).second.snp_reversed = false;

				if(curr_allele_A != first_allele_A || curr_allele_B != first_allele_B) {

					#ifdef DEBUG
						cout << (*it_location).second.rsid << ": alleles "
								<< curr_allele_A << "/" << curr_allele_B << " (Cohort " << i+1 << ")" << " != "
								<< first_allele_A << "/" << first_allele_B << " (Cohort 1)";
					#endif


					if(check_alleles_flipped(first_allele_A, first_allele_B, curr_allele_A, curr_allele_B)) {

						(*it_location).second.snp_flipped = true;
						snp_flipped++;

					}else if(check_alleles_reversed(first_allele_A, first_allele_B, curr_allele_A, curr_allele_B)) {

						(*it_location).second.snp_reversed = true;
						snp_reversed++;

					}else{

						union_pos.erase(it_union);
						snp_removed++;
						break;

					}

				}

			} // if first

			first++;

		} // if(it_location)

	}  // for the number_cohort
}

/* inverse-variance meta analysis */
void META::inverse_variance_method()
{

	/* combine beta and standard error SNP by SNP */
	for(Union::iterator it_union = union_pos.begin(); it_union != union_pos.end(); it_union++)
	{

		double beta_hat = 0.0;
		double se_hat = 0.0;

		combine_beta_se(*it_union, beta_hat, se_hat);
		heterogeneity(meta_snp, beta_hat);

		/* random effect model, variance between group (nu) is taken into account */
		if(methods == 2) combine_beta_se_with_nu(*it_union, beta_hat, se_hat);

		meta_snp.beta = beta_hat;
		meta_snp.se = se_hat;

		#ifdef DEBUG
			cout << "combined beta = " << meta_snp.beta << " combined se = " << meta_snp.se << endl;
		#endif

		/* proper test statistics and p-value */
		double stats =  (beta_hat / se_hat) * (beta_hat / se_hat);
		meta_snp.pvalue = chdtrc(1.0, stats);
		//cout << "stats = " << stats << " pvalue = " << meta_snp.pvalue << endl;

		meta.insert(pair<Positive, MetaSNP>(meta_snp.pos, meta_snp));

		#ifdef DEBUG
			cout << meta_snp.chr << " ";
			cout << meta_snp.rsid << " ";
			cout << meta_snp.pos << " ";
			cout << meta_snp.allele_A << " ";
			cout << meta_snp.allele_B << " ";
			cout << meta_snp.pvalue << " ";
			cout << meta_snp.coded_af << " ";
			cout << meta_snp.beta << " ";
			cout << meta_snp.se << " ";
			cout << meta_snp.Q << " ";
			cout << meta_snp.p_het << " ";
			cout << meta_snp.I2 << endl;
		#endif
	}

}

void META::combine_beta_se(const Positive pos, double& beta_hat, double& se_hat)
{

	beta_i.clear();
	wei_i.clear();

	beta_i.resize(number_cohort);
	wei_i.resize(number_cohort);

	Cohort::iterator it_location;
	pvalue_list.reserve(union_pos.size());

	vector<double> pvalue_i(number_cohort, -1);

	double beta_sum = 0.0;
	double wei_sum = 0.0;
	double coded_count = 0.0;
	double total_count = 0.0;

	int first = 0;
	for(int i = 0; i < number_cohort; i++) {

		it_location = cohort_list[i].find(pos);

		if(it_location != cohort_list[i].end()) {

			if(first == 0) {

				/* rsid, alleles and allele frequenecies */
				if(chr_included) meta_snp.chr = (*it_location).second.chr;
				meta_snp.rsid = (*it_location).second.rsid;
				meta_snp.allele_A = (*it_location).second.allele_A;
				meta_snp.allele_B = (*it_location).second.allele_B;
				meta_snp.pos = (*it_location).first;

			}


			if((*it_location).second.snp_reversed)
			{

				beta_i[i] = -(*it_location).second.beta;
				coded_count += (*it_location).second.genotype_aa + (*it_location).second.genotype_ab / 2.0;

			}else{

				beta_i[i] = (*it_location).second.beta;
				coded_count += (*it_location).second.genotype_bb + (*it_location).second.genotype_ab / 2.0;

			}

			total_count += (*it_location).second.genotype_aa + (*it_location).second.genotype_ab + (*it_location).second.genotype_bb;

			if(gc)
			{

				beta_sum += beta_i[i] / (((*it_location).second.se * (*it_location).second.se) * genomic_controls[i]);
				wei_i[i] = 1.0 / ((*it_location).second.se * (*it_location).second.se * genomic_controls[i]);
				wei_sum += 1.0 / ((*it_location).second.se * (*it_location).second.se * genomic_controls[i]);

			}else{

				beta_sum += beta_i[i] / ((*it_location).second.se * (*it_location).second.se);
				wei_i[i] = 1.0 / ((*it_location).second.se * (*it_location).second.se);
				wei_sum += 1.0 / ((*it_location).second.se * (*it_location).second.se);

			}

			meta_snp.coded_af = coded_count / total_count;
			pvalue_i[i] = (*it_location).second.pvalue;

			first++;

			#ifdef DEBUG
				cout << "cohort " << i + 1 << ": beta = " << beta_i[i] << " se = " << (*it_location).second.se << endl;
			#endif

		}

	}

	pvalue_list.push_back(pvalue_i);

	beta_hat = beta_sum / wei_sum;
	se_hat = sqrt(1.0 / wei_sum);

}


void META::combine_beta_se_with_nu(const Positive pos, double& beta_hat, double& se_hat)
{

	#ifdef DEBUG
		cout << "Q = " << meta_snp.Q << " nu = " << meta_snp.nu << endl;
	#endif

	Cohort::iterator it_location;

	double beta_sum = 0.0;
	double wei_sum = 0.0;

	for(int i = 0; i < number_cohort; i++) {

		it_location = cohort_list[i].find(pos);

		if(it_location != cohort_list[i].end()) {

			if((*it_location).second.snp_reversed) beta_i[i] = -(*it_location).second.beta;
			else beta_i[i] = (*it_location).second.beta;

			if(gc) {

				beta_sum += beta_i[i] / ((*it_location).second.se * (*it_location).second.se * genomic_controls[i] + meta_snp.nu);
				wei_sum += 1.0 / ((*it_location).second.se * (*it_location).second.se * genomic_controls[i] + meta_snp.nu);

			}else{

				beta_sum += beta_i[i] / ((*it_location).second.se * (*it_location).second.se + meta_snp.nu);
				wei_sum += 1.0 / ((*it_location).second.se * (*it_location).second.se + meta_snp.nu);

			}

		}

	}

	beta_hat = beta_sum / wei_sum;
	se_hat = sqrt(1.0 / wei_sum);

}

// meta analysis by combining z-statistics using sample sizes as weights
void META::combine_z_score()
{

  pvalue_list.reserve(union_pos.size());
  for(Union::iterator it_union = union_pos.begin(); it_union != union_pos.end(); it_union++) {

    Cohort::iterator it_location;
    vector<double> pvalue_i(number_cohort, -1.0);
    vector<double> wei_i(number_cohort, 0.0);
    vector<double> beta_i(number_cohort, 0.0);
    vector<double> se_i(number_cohort, 0.0);
    vector<double> stat_i(number_cohort, 0.0);
    double wei_sum = 0.0;
    double coded_count = 0.0;
    double total_count = 0.0;

    for(int i = 0; i < number_cohort; i++)
    {

    	it_location = cohort_list[i].find(*it_union);

    	if(it_location != cohort_list[i].end())
    	{

    		/* rsid, alleles and allele frequenecies */
    		meta_snp.rsid = (*it_location).second.rsid;
    		meta_snp.allele_A = (*it_location).second.allele_A;
    		meta_snp.allele_B = (*it_location).second.allele_B;
    		meta_snp.pos = (*it_location).first;

    		if((*it_location).second.snp_reversed) beta_i[i] = -(*it_location).second.beta;
    		else beta_i[i] = (*it_location).second.beta;

    		se_i[i] = (*it_location).second.se;

    		if(gc) stat_i[i] = beta_i[i] / (se_i[i] * sqrt(genomic_controls[i]));
    		else stat_i[i] = beta_i[i] / se_i[i];

    		coded_count += (*it_location).second.genotype_aa + (*it_location).second.genotype_ab / 2.0;

    		total_count += (*it_location).second.genotype_aa + (*it_location).second.genotype_ab + (*it_location).second.genotype_bb;

    		wei_i[i] = sample_sizes[i];
    		wei_sum += wei_i[i];

    		#ifdef DEBUG
			cout << "Cohort " << i + 1 << ": beta = " << beta_i[i] << " se = " << se_i[i] << " wei = " << wei_i[i] << endl;
		#endif

    		pvalue_i[i] = (*it_location).second.pvalue;

    	}

	}

    double stats_meta = 0.0;
    for(int i = 0; i < number_cohort; i++) stats_meta += stat_i[i] * sqrt(wei_i[i] / wei_sum);

    pvalue_list.push_back(pvalue_i);
    meta_snp.pvalue = 2.0 * ndtr(-abs(stats_meta));

    #ifdef DEBUG
	cout << "Statistics = " << stats_meta << " P-value = " << meta_snp.pvalue << endl;
    #endif

    meta.insert(pair<Positive, MetaSNP>(meta_snp.pos, meta_snp));

  }  // for union

}


bool META::check_alleles_flipped(char first_allele_A, char first_allele_B, char curr_allele_A, char curr_allele_B)
{

	const char alleles[] = "ATGC";

	bool flip = false;
	if(curr_allele_A != first_allele_B && curr_allele_B != first_allele_A)
	{

		bool flip1 = false;
		for(int k = 0; k < 4; k++) {

			if(curr_allele_A == alleles[k]) {

				if(k < 2) {

					if(first_allele_A == alleles[1 - k]) flip1 = true;

				}else{

					if(first_allele_A == alleles[5 - k]) flip1 = true;

				}

			}

		}

		bool flip2 = false;
		for(int k = 0; k < 4; k++) {

			if(curr_allele_B == alleles[k]) {

				if(k < 2) {

					if(first_allele_B == alleles[1 - k]) flip2 = true;

				}else{

					if(first_allele_B == alleles[5 - k]) flip2 = true;

				}

			}

		}

		flip = flip1 * flip2;
	}

	return flip;

}


bool META::check_alleles_reversed(char first_allele_A, char first_allele_B, char curr_allele_A, char curr_allele_B)
{

	bool reverse = false;

	if(first_allele_A == curr_allele_B && first_allele_B == curr_allele_A) reverse = true;

	else{

		const char alleles[] = "ATGC";

		bool reverse1 = false;

		for(int k = 0; k < 4; k++) {

			if(curr_allele_A == alleles[k]) {

				if(k < 2) {

					if(first_allele_B == alleles[1 - k]) reverse1 = true;

				}else{

					if(first_allele_B == alleles[5 - k]) reverse1 = true;

				}

			}

		}

		bool reverse2 = false;

		for(int k = 0; k < 4; k++) {

			if(curr_allele_B == alleles[k]) {

				if(k < 2) {

					if(first_allele_A == alleles[1 - k]) reverse2 = true;

				}else{

					if(first_allele_A == alleles[5 - k]) reverse2 = true;

				}

			}

		}

		reverse = reverse1 * reverse2;

	}

	return reverse;

}


void META::heterogeneity(MetaSNP& meta_snp, const double beta_hat)
{

	meta_snp.Q = 0;
	int degree = 0;
	double beta_sum = 0.0;
	double wei_sum = 0.0;
	double wei2_sum = 0.0;

	for(int i = 0; i < number_cohort; i++) {

		meta_snp.Q += (wei_i[i] * (beta_i[i] - beta_hat) * (beta_i[i] - beta_hat));
		beta_sum += beta_i[i];
		wei_sum += wei_i[i];
		wei2_sum += wei_i[i] * wei_i[i];
		if(wei_i[i] != 0) degree++;

	}

	//cout << "degree = " << degree << " wei2 = " << wei2_sum << "wei =  " << wei_sum<< endl;

	if(degree > 1) {

		if(meta_snp.Q > degree - 1) meta_snp.nu = (meta_snp.Q - degree + 1) / (wei_sum - wei2_sum / wei_sum);
		else meta_snp.nu = 0.0;
		//cout << "id = " << meta_snp.rsid << " beta.hat = " << beta_hat << " Q = " << meta_snp.Q << endl;
		meta_snp.p_het = chdtrc(double (degree - 1), meta_snp.Q);
		meta_snp.I2 = 100 * double (meta_snp.Q - degree + 1) / meta_snp.Q;
		if(meta_snp.I2 < 0) meta_snp.I2 = 0;

	}else{

		meta_snp.nu = 0;
		meta_snp.p_het = -1;
		meta_snp.I2 = -1;

	}

}
