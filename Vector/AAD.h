#pragma once
#pragma once
#ifndef AADH
#define AADH
#include <cmath>
#include <vector>
#include <autodiff.h>

using namespace std;

class aad
{
private:

	class adjoint
	{
	private:
		friend class aad;

		double  value;
		vector<aad> tape;
		vector<double> weights;
		double grad;
	public:
		adjoint();
		adjoint(double);
		adjoint(const adjoint&);

		double GetValue();
		void  SetValue(double v);

		double GetGrad();
		void SetGrad(double g);

		vector<aad> GetTape();
		void SetTape(vector<aad> t);
		void AppendToTape(aad& a);

		vector<double> GetWeights();
		void SetWeights(vector<double> w);
		void AppendToWeights(double a);

		void adjoint::operator = (adjoint &c);
	};

	adjoint * var;
	
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

	friend aad operator * (const double c, aad &f);
	friend aad operator / (const double c, aad &f);

	friend aad operator + (aad &c, aad &f);
	friend aad operator + (const double c, aad &f);
	friend aad operator + (aad &f, const double c);
	friend aad operator + (aad &c);

	friend aad operator - (aad &c, aad &f);
	friend aad operator - (const double c, aad &f);
	friend aad operator - (aad &f, const double c);
	friend aad operator - (aad &c);



	friend aad sin(aad &f);
	friend aad cos(aad &f);
	friend aad log(aad &f);
	friend aad exp(aad &f);
	friend aad sqrt(aad &f);
	friend aad N(aad &f); // standart normal distribution
	friend double N(const double &f); // standart normal distribution 
	
public:
	aad();
	aad(double);
	aad(const adjoint&);
	~aad();
	void operator = (aad &c);
	void operator = (double);

	aad operator += (aad &);
	aad operator += (const double&);

	aad operator * (aad &);
	aad  operator * (const double&);

	aad operator / (aad &);
	aad operator / (const double&);

	bool operator > (const aad &) const;
	bool operator >= (const aad &) const;
	bool operator < (const aad &) const;
	bool operator <= (const aad &) const;

	bool operator > (const double &) const;
	bool operator >= (const double &) const;
	bool operator < (const double &) const;
	bool operator <= (const double &) const;

	adjoint* GetVar();

	double GetGradient();
	void SetGradient(double g);

	void Reset();
};

#endif

