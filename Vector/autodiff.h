#define _USE_MATH_DEFINES
#pragma once
#ifndef autodiffH
#define autodiffH
#include <math.h>
#include <vector>

using namespace std;

class autodiff
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

	friend autodiff operator * (const double c, const autodiff &f);
	friend autodiff operator / (const double c, const autodiff &f);

	friend autodiff operator + (const autodiff &c, const autodiff &f);
	friend autodiff operator + (const double c, const autodiff &f);
	friend autodiff operator + (const autodiff &f, const double c);
	friend autodiff operator + (const autodiff &c); 

	friend autodiff operator - (const autodiff &c, const autodiff &f);
	friend autodiff operator - (const double c, const autodiff &f);
	friend autodiff operator - (const autodiff &f, const double c);
	friend autodiff operator - (const autodiff &c); 
	
	friend autodiff sin(autodiff &f);
	friend autodiff cos(autodiff &f);
	friend autodiff log(autodiff &f);
	friend autodiff exp(autodiff &f);
	friend autodiff sqrt(autodiff &f);
	friend autodiff N(autodiff &f); // standart normal distribution
	friend double N(const double &f); // standart normal distribution 
	friend double dN(const double &f); // density

public:
	double  value; 
	vector<double>  deriv;
	autodiff(double, vector<double>); 
	autodiff();
	autodiff(double, int);
	autodiff(double);
	void operator = (autodiff &c);
	void operator = (double);
	
	autodiff operator * (const autodiff &);
	autodiff  operator * (const double &);

	autodiff operator / (const autodiff &);
	autodiff operator / (const double &);

	bool operator > (const autodiff &) const;
	bool operator >= (const autodiff &) const;
	bool operator < (const autodiff &) const;
	bool operator <= (const autodiff &) const;

	bool operator > (const double &) const;
	bool operator >= (const double &) const;
	bool operator < (const double &) const;
	bool operator <= (const double &) const;
};

#endif

