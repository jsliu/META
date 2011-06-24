/*
 * error.h
 *
 *  Created on: 13 Jul 2010
 *      Author: jsliu
 */

#ifndef ERROR_H_
#define ERROR_H_

#include<stdexcept>

using namespace std;

class BadFile: public exception
{
public:
	BadFile(const string error): msg(error.c_str()) {};
	// the standard exception has a virtual destructor promised
	// not to throw any exception, so the derived class shouldn't
	// throw any exception; the keyword virtual is not necessary
	virtual ~BadFile() throw() {};
	// return the error message without throwing any exception
	virtual const char* what() const throw() { return msg; };

private:
	const char* msg;
};

#endif /* ERROR_H_ */
