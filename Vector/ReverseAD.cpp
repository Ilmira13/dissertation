#define _USE_MATH_DEFINES
#pragma hdrstop
#include <ReverseAD.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

// ------------------------ adjoint --------------------------------------------------------------------------
ReverseAD::adjoint::adjoint(double a)
{
	value = a;
	tape.clear();
	weights.clear();
	grad = NAN;
}

ReverseAD::adjoint::adjoint()
{
	value = 0;
	tape.clear();
	weights.clear();
	grad = NAN;
}

ReverseAD::adjoint::adjoint(const adjoint &c)
{
	value = c.value;
	tape = c.tape;
	weights = c.weights;
	grad = c.grad;
}

void ReverseAD::adjoint::operator = (adjoint &c)
{
	this->value = c.value;
	this->tape = c.tape;
	this->weights = c.weights;
	this->grad = c.grad;
}

double ReverseAD::adjoint::GetValue()
{
	return value;
}

void  ReverseAD::adjoint::SetValue(double v)
{
	value = v;
}

double ReverseAD::adjoint::GetGrad()
{
	return grad;
}

void ReverseAD::adjoint::SetGrad(double g)
{
	grad = g;
}

vector<ReverseAD> ReverseAD::adjoint::GetTape()
{
	return tape;
}

void ReverseAD::adjoint::SetTape(vector<ReverseAD> t)
{
	tape = t;
}

void ReverseAD::adjoint::AppendToTape(ReverseAD& a)
{
	tape.push_back(a);
}

vector<double> ReverseAD::adjoint::GetWeights()
{
	return weights;
}

void ReverseAD::adjoint::SetWeights(vector<double> w)
{
	weights = w;
}

void ReverseAD::adjoint::AppendToWeights(double a)
{
	weights.push_back(a);
}

// --------------------------------------------------------------------------------------------------

// ---------------- ReverseAD ----------------------------------------------------------------------------------
ReverseAD::ReverseAD()
{
	var = new adjoint(0.0);
}

ReverseAD::ReverseAD(double a)
{
	var = new adjoint(a);
}

ReverseAD::ReverseAD(const adjoint &a)
{
	var = new adjoint(a);
}

/*ReverseAD::ReverseAD(const ReverseAD &a)
{
	var = new adjoint(*a.var);
}

ReverseAD::ReverseAD(ReverseAD &a)
{
	var = new adjoint(*a.var);
}*/

ReverseAD::~ReverseAD()
{
	//delete var;
}


void ReverseAD::operator = (ReverseAD &c)
{
	//delete this->var;
	this->var = new adjoint(c.var->value);
	c.var->weights.push_back(1.0);
	c.var->tape.push_back(*this);
	//this->var = c.var;
}

void ReverseAD::operator=(double c)
{
	//delete this->var;
	this->var = new adjoint(c);
}


ReverseAD ReverseAD::operator += (ReverseAD &a)
{
	//ReverseAD res(this->var->value + a.var->value);
	//this->var->weights.push_back(1.0);
	//this->var->tape.push_back(res);
	//a.var->weights.push_back(1.0);
	//a.var->tape.push_back(res);
	this->var->value += a.var->value;
	a.var->weights.push_back(1.0);
	a.var->tape.push_back(*this);
	return *this;
}

ReverseAD ReverseAD::operator += (const double &c)
{
	this->var->value += c;
	return *this;
}


ReverseAD ReverseAD::operator *= (ReverseAD &a)
{
	//ReverseAD res(this->var->value + a.var->value);
	//this->var->weights.push_back(1.0);
	//this->var->tape.push_back(res);
	//a.var->weights.push_back(1.0);
	//a.var->tape.push_back(res);
	this->var->value *= a.var->value;
	a.var->weights.push_back(a.var->value);
	a.var->tape.push_back(*this);
	return *this;
}

ReverseAD ReverseAD::operator *= (const double &c)
{
	this->var->value *= c;
	this->var->weights.push_back(c);
	this->var->tape.push_back(*this);
	return *this;
}


ReverseAD operator + (ReverseAD  &c)
{
	return c;
}

ReverseAD  operator - (ReverseAD  &c)
{
	ReverseAD res(-c.var->GetValue());
	res.var->SetWeights(-c.var->GetWeights());
	if (!isnan(res.var->GetGrad()))
		res.var->SetGrad(-res.var->GetGrad());
	return res;
}

ReverseAD  ReverseAD::operator * (ReverseAD &c)
{
	ReverseAD res(this->var->value * c.var->value);
	this->var->weights.push_back(c.var->value);
	this->var->tape.push_back(res);
	c.var->weights.push_back(this->var->value);
	c.var->tape.push_back(res);
	return res;
}
ReverseAD ReverseAD::operator * (const double &c)
{
	ReverseAD res(this->var->value * c);
	this->var->weights.push_back(c);
	this->var->tape.push_back(res);
	return res;
}

ReverseAD operator * (const double c, ReverseAD &f)
{
	ReverseAD res(f.var->GetValue() * c);
	f.var->AppendToWeights(c);
	f.var->AppendToTape(res);
	return res;
}
/*ReverseAD operator * (ReverseAD &f, const double c)
{
	ReverseAD res(f.var->GetValue() * c);
	f.var->AppendToWeights(c);
	f.var->AppendToTape(res);
	return res;
}
ReverseAD operator * (ReverseAD &a, ReverseAD &f)
{
	ReverseAD res(f.var->GetValue() * a.var->GetValue());
	a.var->AppendToWeights(f.var->GetValue());
	a.var->AppendToTape(res);
	f.var->AppendToWeights(a.var->GetValue());
	f.var->AppendToTape(res);
	return res;
}*/

ReverseAD  ReverseAD::operator / (ReverseAD &c)
{
	ReverseAD res(this->var->value / c.var->value);
	this->var->weights.push_back(1.0 / c.var->value);
	this->var->tape.push_back(res);
	c.var->weights.push_back(- this->var->value / c.var->value / c.var->value);
	c.var->tape.push_back(res);
	return res;
}
ReverseAD ReverseAD::operator / (const double &c)
{
	ReverseAD res(this->var->value / c);
	this->var->weights.push_back(1.0 / c);
	this->var->tape.push_back(res);
	return res;
}

ReverseAD operator / (const double c, ReverseAD &f)
{
	ReverseAD res(c / f.GetVar()->GetValue());
	f.GetVar()->AppendToWeights(-c / f.GetVar()->GetValue() / f.GetVar()->GetValue());
	f.GetVar()->AppendToTape(res);
	return res;
}
/*ReverseAD operator / (ReverseAD &f, const double c)
{
	ReverseAD res(c / f.GetVar()->GetValue());
	f.GetVar()->AppendToWeights(-c / f.GetVar()->GetValue() / f.GetVar()->GetValue());
	f.GetVar()->AppendToTape(res);
	return res;
}
ReverseAD operator / (ReverseAD &a, ReverseAD &f)
{
	ReverseAD res(a.GetVar()->GetValue() / f.GetVar()->GetValue());
	a.GetVar()->AppendToWeights(1.0 / f.GetVar()->GetValue());
	a.GetVar()->AppendToTape(res);
	f.GetVar()->AppendToWeights(-a.GetVar()->GetValue() / f.GetVar()->GetValue() / f.GetVar()->GetValue());
	f.GetVar()->AppendToTape(res);
	return res;
}*/

ReverseAD operator + (ReverseAD &c, ReverseAD &f)
{
	ReverseAD res(f.GetVar()->GetValue() + c.GetVar()->GetValue());
	f.GetVar()->AppendToWeights(1.0);
	f.GetVar()->AppendToTape(res);
	c.GetVar()->AppendToWeights(1.0);
	c.GetVar()->AppendToTape(res);
	return res;
}

ReverseAD operator + (ReverseAD &f, const double c)
{
	ReverseAD res(f.GetVar()->GetValue() + c);
	f.GetVar()->AppendToWeights(1.0);
	f.GetVar()->AppendToTape(res);
	return res;
}

ReverseAD operator + (const double c, ReverseAD &f)
{
	ReverseAD res(f.GetVar()->GetValue() + c);
	f.GetVar()->AppendToWeights(1.0);
	f.GetVar()->AppendToTape(res);
	return res;
}


ReverseAD operator - (ReverseAD &c, ReverseAD &f)
{
	ReverseAD res(c.GetVar()->GetValue() - f.GetVar()->GetValue());
	f.GetVar()->AppendToWeights(-1.0);
	f.GetVar()->AppendToTape(res);
	c.GetVar()->AppendToWeights(1.0);
	c.GetVar()->AppendToTape(res);
	return res;
}
ReverseAD operator - (ReverseAD &f, const double c)
{
	ReverseAD res(f.GetVar()->GetValue() - c);
	f.GetVar()->AppendToWeights(1.0);
	f.GetVar()->AppendToTape(res);
	return res;
}
ReverseAD operator - (const double c, ReverseAD &f)
{
	ReverseAD res(-f.GetVar()->GetValue() + c);
	f.GetVar()->AppendToWeights(-1.0);
	f.GetVar()->AppendToTape(res);
	return res;
}

bool ReverseAD::operator>(const ReverseAD & c) const
{
	return this->var->GetValue() > c.var->GetValue();
}
bool ReverseAD::operator>=(const ReverseAD & c) const
{
	return this->var->GetValue() >= c.var->GetValue();
}
bool ReverseAD::operator<(const ReverseAD & c) const
{
	return this->var->GetValue() < c.var->GetValue();
}
bool ReverseAD::operator<=(const ReverseAD & c) const
{
	return this->var->GetValue() <= c.var->GetValue();
}
bool ReverseAD::operator>(const double & c) const
{
	return this->var->GetValue() > c;
}
bool ReverseAD::operator>=(const double & c) const
{
	return this->var->GetValue() >= c;
}
bool ReverseAD::operator<(const double & c) const
{
	return this->var->GetValue() < c;
}
bool ReverseAD::operator<=(const double & c) const
{
	return this->var->GetValue() <= c;
}

ReverseAD::adjoint* ReverseAD::GetVar()
{
	return var;
}
double ReverseAD::GetGradient()
{
	if (var->weights.size() > 0 && isnan(var->grad))
	{
		var->grad = 0.0;
		for (size_t i = 0; i < var->weights.size(); ++i)
		{
			var->grad += var->weights[i] * var->tape[i].GetGradient();
		}
	}
	else if (var->weights.size() == 0 && isnan(var->grad))
		var->grad = 1.0;
	return var->grad;
}
void ReverseAD::SetGradient(double g)
{
	var->grad = g;
}

void ReverseAD::Reset()
{
	double val = var->value;
	//delete var;
	var = new adjoint(val);
}

ReverseAD N(ReverseAD &f) // standard normal distribution
{
	ReverseAD res(0.5*(1 + erf(f.GetVar()->GetValue() / sqrt(2))));
	f.GetVar()->AppendToWeights(1.0 / sqrt(2.0 * M_PI)*exp(-0.5*pow(f.GetVar()->GetValue(), 2)));
	f.GetVar()->AppendToTape(res);
	return res;
}

ReverseAD sin(ReverseAD &f)
{
	ReverseAD res(sin(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(cos(f.GetVar()->GetValue()));
	f.GetVar()->AppendToTape(res);
	return res;
}

ReverseAD cos(ReverseAD &f)
{
	ReverseAD res(cos(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(-sin(f.GetVar()->GetValue()));
	f.GetVar()->AppendToTape(res);
	return res;
}

ReverseAD log(ReverseAD &f)
{
	ReverseAD res(log(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(1.0 / f.GetVar()->GetValue());
	f.GetVar()->AppendToTape(res);
	return res;
}

ReverseAD exp(ReverseAD &f)
{
	ReverseAD res(exp(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(exp(f.GetVar()->GetValue()));
	f.GetVar()->AppendToTape(res);
	return res;
}

ReverseAD sqrt(ReverseAD &f)
{
	ReverseAD res(sqrt(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(-0.5 / sqrt(f.GetVar()->GetValue()));
	f.GetVar()->AppendToTape(res);
	return res;
}







