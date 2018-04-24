#pragma once
#pragma once
#ifndef ReverseADH
#define ReverseADH
#include <cmath>
#include <vector>
#include <ForwardAD.h>

using namespace std;

class ReverseAD
{
public:

	class adjoint
	{
	public:
		friend class ReverseAD;

		double  value;
		vector<ReverseAD> tape;
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

		vector<ReverseAD> GetTape();
		void SetTape(vector<ReverseAD> t);
		void AppendToTape(ReverseAD& a);

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

	friend ReverseAD operator * (const double c, ReverseAD &f);
	//friend ReverseAD operator * (ReverseAD &f, const double c);
	//friend ReverseAD operator * (ReverseAD &a, ReverseAD &f);
	friend ReverseAD operator / (const double c, ReverseAD &f);
	//friend ReverseAD operator / (ReverseAD &f, const double c);
	//friend ReverseAD operator / (ReverseAD &a, ReverseAD &f);

	friend ReverseAD operator + (ReverseAD &c, ReverseAD &f);
	friend ReverseAD operator + (const double c, ReverseAD &f);
	friend ReverseAD operator + (ReverseAD &f, const double c);
	friend ReverseAD operator + (ReverseAD &c);

	friend ReverseAD operator - (ReverseAD &c, ReverseAD &f);
	friend ReverseAD operator - (const double c, ReverseAD &f);
	friend ReverseAD operator - (ReverseAD &f, const double c);
	friend ReverseAD operator - (ReverseAD &c);



	friend ReverseAD sin(ReverseAD &f);
	friend ReverseAD cos(ReverseAD &f);
	friend ReverseAD log(ReverseAD &f);
	friend ReverseAD exp(ReverseAD &f);
	friend ReverseAD sqrt(ReverseAD &f);
	friend ReverseAD N(ReverseAD &f); // standart normal distribution
	friend double N(const double &f); // standart normal distribution 
	
public:
	ReverseAD();
	ReverseAD(double);
	ReverseAD(const adjoint&);
	//ReverseAD(const ReverseAD &a);
	//ReverseAD(ReverseAD &a);
	~ReverseAD();
	void operator = (ReverseAD &c);
	void operator = (double);

	ReverseAD operator += (ReverseAD &); // does not work
	ReverseAD operator += (const double&); // does not work

	ReverseAD operator *= (ReverseAD &); // does not work
	ReverseAD operator *= (const double&); // does not work

	ReverseAD operator * (ReverseAD &);
	ReverseAD  operator * (const double&);

	ReverseAD operator / (ReverseAD &);
	ReverseAD operator / (const double&);

	bool operator > (const ReverseAD &) const;
	bool operator >= (const ReverseAD &) const;
	bool operator < (const ReverseAD &) const;
	bool operator <= (const ReverseAD &) const;

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

