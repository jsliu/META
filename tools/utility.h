#ifndef UTILITY_H_
#define UTILITY_H_

#include <string>
#include <limits>
#include <sys/stat.h>

using namespace std;

template<class T> inline bool is_nan(T value)
{
	return value != value;
}

/* to see if the value is infinite */
template<class T> inline bool is_inf(T value)
{
	return value >= std::numeric_limits<T>::min() && value <= std::numeric_limits<T>::max();
}


inline bool file_exists(string strFilename)
{
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0)
  {
	  // We were able to get the file attributes
      // so the file obviously exists.
      blnReturn = true;

  } else {

	  // We were not able to get the file attributes.
      // This may mean that we don't have permission to
      // access the folder which contains this file. If you
      // need to do that level of checking, lookup the
      // return values of stat which will give you
      // more details on why stat failed.
      blnReturn = false;

  }

  return(blnReturn);
}


inline char string2char(const string& s)
{
	istringstream iss(s, istringstream::in);
	char c = iss.get();

	return c;

}


inline string char2string(const char* ch)
{
	string str(ch);
	return str;
}


inline int string2int(const string& s)
{
	istringstream iss(s, istringstream::in);
	int i;
	iss >> i;

	return i;
}

inline string int2string(const int i)
{
	ostringstream oss;
	oss << i;
	string s = oss.str();

	return s;
}

inline long string2long(const string& s)
{
	istringstream iss(s, istringstream::in);
	long i;
	iss >> i;

	return i;
}


inline double string2double(const string& s)
{

	istringstream iss(s, istringstream::in);
	double d;
	iss >> d;

	return d;

}


inline vector<double> string2double(const vector<string>& s)
{
	vector<double> d;
	d.reserve(s.size());

	for(size_t i = 0; i < s.size(); i++)
	{
		istringstream iss(s[i], istringstream::in);
		double tmp;
		iss >> tmp;
		d.push_back(tmp);
	}

	return d;
}


inline string double2string(const double d)
{

	ostringstream oss;
	oss << d;
	string s = oss.str();

	return s;
}

#endif
