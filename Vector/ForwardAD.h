#define _USE_MATH_DEFINES
#pragma once
#ifndef ForwardADH
#define ForwardADH
#include <math.h>
#include <vector>

using namespace std;

class ForwardAD
{
	friend vector<double> operator + (const vector<double>u); 
	friend vector<double> operator - (const vector<double>u);
    friend vector<double> operator + (const vector<double>u, const vector<double> h);
	friend vector<double> operator - (const vector<double>u, const vector<double> h);
	friend vector<double> operator * (const vector<double>u, const vector<double> h);
	friend vector<double> operator * (const double &u, const vector<double> h);
	friend vector<double> operator * (const vector<double> h, const double &u);
	friend vector<double> operator / (const vector<double>u, const vector<double> h);
	friend vector<double> operator / (const double &u, const vector<double> h);
	friend vector<double> operator / (const vector<double> h, const double &u);

	friend ForwardAD operator * (const double c, const ForwardAD &f);
	friend ForwardAD operator / (const double c, const ForwardAD &f);

	friend ForwardAD operator + (const ForwardAD &c, const ForwardAD &f);
	friend ForwardAD operator + (const double c, const ForwardAD &f);
	friend ForwardAD operator + (const ForwardAD &f, const double c);
	friend ForwardAD operator + (const ForwardAD &c); 

	friend ForwardAD operator - (const ForwardAD &c, const ForwardAD &f);
	friend ForwardAD operator - (const double c, const ForwardAD &f);
	friend ForwardAD operator - (const ForwardAD &f, const double c);
	friend ForwardAD operator - (const ForwardAD &c); 
	
	friend ForwardAD sin(ForwardAD &f);
	friend ForwardAD cos(ForwardAD &f);
	friend ForwardAD log(ForwardAD &f);
	friend ForwardAD exp(ForwardAD &f);
	friend ForwardAD sqrt(ForwardAD &f);
	friend ForwardAD N(ForwardAD &f); // standart normal distribution
	friend double N(const double &f); // standart normal distribution 
	friend double dN(const double &f); // density

public:
	double  value; 
	vector<double>  deriv;
	ForwardAD(double, vector<double>); 
	ForwardAD();
	ForwardAD(double, int);
	ForwardAD(double);
	void operator = (ForwardAD &c);
	void operator = (double);
	
	ForwardAD operator * (const ForwardAD &);
	ForwardAD  operator * (const double &);

	ForwardAD operator / (const ForwardAD &);
	ForwardAD operator / (const double &);

	bool operator > (const ForwardAD &) const;
	bool operator >= (const ForwardAD &) const;
	bool operator < (const ForwardAD &) const;
	bool operator <= (const ForwardAD &) const;

	bool operator > (const double &) const;
	bool operator >= (const double &) const;
	bool operator < (const double &) const;
	bool operator <= (const double &) const;
};

#endif

