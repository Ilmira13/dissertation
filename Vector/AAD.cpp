#define _USE_MATH_DEFINES
#pragma hdrstop
#include <aad.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

// ------------------------ adjoint --------------------------------------------------------------------------
aad::adjoint::adjoint(double a)
{
	value = a;
	tape.clear();
	weights.clear();
	grad = NAN;
}

aad::adjoint::adjoint()
{
	value = 0;
	tape.clear();
	weights.clear();
	grad = NAN;
}

aad::adjoint::adjoint(const adjoint &c)
{
	value = c.value;
	tape = c.tape;
	weights = c.weights;
	grad = c.grad;
}

void aad::adjoint::operator = (adjoint &c)
{
	this->value = c.value;
	this->tape = c.tape;
	this->weights = c.weights;
	this->grad = c.grad;
}

double aad::adjoint::GetValue()
{
	return value;
}

void  aad::adjoint::SetValue(double v)
{
	value = v;
}

double aad::adjoint::GetGrad()
{
	return grad;
}

void aad::adjoint::SetGrad(double g)
{
	grad = g;
}

vector<aad> aad::adjoint::GetTape()
{
	return tape;
}

void aad::adjoint::SetTape(vector<aad> t)
{
	tape = t;
}

void aad::adjoint::AppendToTape(aad& a)
{
	tape.push_back(a);
}

vector<double> aad::adjoint::GetWeights()
{
	return weights;
}

void aad::adjoint::SetWeights(vector<double> w)
{
	weights = w;
}

void aad::adjoint::AppendToWeights(double a)
{
	weights.push_back(a);
}

// --------------------------------------------------------------------------------------------------

// ---------------- aad ----------------------------------------------------------------------------------
aad::aad()
{
	var = new adjoint(0.0);
}

aad::aad(double a)
{
	var = new adjoint(a);
}

aad::aad(const adjoint &a)
{
	var = new adjoint(a);
}

/*aad::aad(const aad &a)
{
	var = new adjoint(*a.var);
}

aad::aad(aad &a)
{
	var = new adjoint(*a.var);
}*/

aad::~aad()
{
	//delete var;
}


void aad::operator = (aad &c)
{
	//delete this->var;
	this->var = new adjoint(c.var->value);
	c.var->weights.push_back(1.0);
	c.var->tape.push_back(*this);
	//this->var = c.var;
}

void aad::operator=(double c)
{
	//delete this->var;
	this->var = new adjoint(c);
}


aad aad::operator += (aad &a)
{
	//aad res(this->var->value + a.var->value);
	//this->var->weights.push_back(1.0);
	//this->var->tape.push_back(res);
	//a.var->weights.push_back(1.0);
	//a.var->tape.push_back(res);
	this->var->value += a.var->value;
	a.var->weights.push_back(1.0);
	a.var->tape.push_back(*this);
	return *this;
}

aad aad::operator += (const double &c)
{
	this->var->value += c;
	return *this;
}


aad aad::operator *= (aad &a)
{
	//aad res(this->var->value + a.var->value);
	//this->var->weights.push_back(1.0);
	//this->var->tape.push_back(res);
	//a.var->weights.push_back(1.0);
	//a.var->tape.push_back(res);
	this->var->value *= a.var->value;
	a.var->weights.push_back(a.var->value);
	a.var->tape.push_back(*this);
	return *this;
}

aad aad::operator *= (const double &c)
{
	this->var->value *= c;
	this->var->weights.push_back(c);
	this->var->tape.push_back(*this);
	return *this;
}


aad operator + (aad  &c)
{
	return c;
}

aad  operator - (aad  &c)
{
	aad res(-c.var->GetValue());
	res.var->SetWeights(-c.var->GetWeights());
	if (!isnan(res.var->GetGrad()))
		res.var->SetGrad(-res.var->GetGrad());
	return res;
}

aad  aad::operator * (aad &c)
{
	aad res(this->var->value * c.var->value);
	this->var->weights.push_back(c.var->value);
	this->var->tape.push_back(res);
	c.var->weights.push_back(this->var->value);
	c.var->tape.push_back(res);
	return res;
}
aad aad::operator * (const double &c)
{
	aad res(this->var->value * c);
	this->var->weights.push_back(c);
	this->var->tape.push_back(res);
	return res;
}

aad operator * (const double c, aad &f)
{
	aad res(f.var->GetValue() * c);
	f.var->AppendToWeights(c);
	f.var->AppendToTape(res);
	return res;
}
/*aad operator * (aad &f, const double c)
{
	aad res(f.var->GetValue() * c);
	f.var->AppendToWeights(c);
	f.var->AppendToTape(res);
	return res;
}
aad operator * (aad &a, aad &f)
{
	aad res(f.var->GetValue() * a.var->GetValue());
	a.var->AppendToWeights(f.var->GetValue());
	a.var->AppendToTape(res);
	f.var->AppendToWeights(a.var->GetValue());
	f.var->AppendToTape(res);
	return res;
}*/

aad  aad::operator / (aad &c)
{
	aad res(this->var->value / c.var->value);
	this->var->weights.push_back(1.0 / c.var->value);
	this->var->tape.push_back(res);
	c.var->weights.push_back(- this->var->value / c.var->value / c.var->value);
	c.var->tape.push_back(res);
	return res;
}
aad aad::operator / (const double &c)
{
	aad res(this->var->value / c);
	this->var->weights.push_back(1.0 / c);
	this->var->tape.push_back(res);
	return res;
}

aad operator / (const double c, aad &f)
{
	aad res(c / f.GetVar()->GetValue());
	f.GetVar()->AppendToWeights(-c / f.GetVar()->GetValue() / f.GetVar()->GetValue());
	f.GetVar()->AppendToTape(res);
	return res;
}
/*aad operator / (aad &f, const double c)
{
	aad res(c / f.GetVar()->GetValue());
	f.GetVar()->AppendToWeights(-c / f.GetVar()->GetValue() / f.GetVar()->GetValue());
	f.GetVar()->AppendToTape(res);
	return res;
}
aad operator / (aad &a, aad &f)
{
	aad res(a.GetVar()->GetValue() / f.GetVar()->GetValue());
	a.GetVar()->AppendToWeights(1.0 / f.GetVar()->GetValue());
	a.GetVar()->AppendToTape(res);
	f.GetVar()->AppendToWeights(-a.GetVar()->GetValue() / f.GetVar()->GetValue() / f.GetVar()->GetValue());
	f.GetVar()->AppendToTape(res);
	return res;
}*/

aad operator + (aad &c, aad &f)
{
	aad res(f.GetVar()->GetValue() + c.GetVar()->GetValue());
	f.GetVar()->AppendToWeights(1.0);
	f.GetVar()->AppendToTape(res);
	c.GetVar()->AppendToWeights(1.0);
	c.GetVar()->AppendToTape(res);
	return res;
}

aad operator + (aad &f, const double c)
{
	aad res(f.GetVar()->GetValue() + c);
	f.GetVar()->AppendToWeights(1.0);
	f.GetVar()->AppendToTape(res);
	return res;
}

aad operator + (const double c, aad &f)
{
	aad res(f.GetVar()->GetValue() + c);
	f.GetVar()->AppendToWeights(1.0);
	f.GetVar()->AppendToTape(res);
	return res;
}


aad operator - (aad &c, aad &f)
{
	aad res(c.GetVar()->GetValue() - f.GetVar()->GetValue());
	f.GetVar()->AppendToWeights(-1.0);
	f.GetVar()->AppendToTape(res);
	c.GetVar()->AppendToWeights(1.0);
	c.GetVar()->AppendToTape(res);
	return res;
}
aad operator - (aad &f, const double c)
{
	aad res(f.GetVar()->GetValue() - c);
	f.GetVar()->AppendToWeights(1.0);
	f.GetVar()->AppendToTape(res);
	return res;
}
aad operator - (const double c, aad &f)
{
	aad res(-f.GetVar()->GetValue() + c);
	f.GetVar()->AppendToWeights(-1.0);
	f.GetVar()->AppendToTape(res);
	return res;
}

bool aad::operator>(const aad & c) const
{
	return this->var->GetValue() > c.var->GetValue();
}
bool aad::operator>=(const aad & c) const
{
	return this->var->GetValue() >= c.var->GetValue();
}
bool aad::operator<(const aad & c) const
{
	return this->var->GetValue() < c.var->GetValue();
}
bool aad::operator<=(const aad & c) const
{
	return this->var->GetValue() <= c.var->GetValue();
}
bool aad::operator>(const double & c) const
{
	return this->var->GetValue() > c;
}
bool aad::operator>=(const double & c) const
{
	return this->var->GetValue() >= c;
}
bool aad::operator<(const double & c) const
{
	return this->var->GetValue() < c;
}
bool aad::operator<=(const double & c) const
{
	return this->var->GetValue() <= c;
}

aad::adjoint* aad::GetVar()
{
	return var;
}
double aad::GetGradient()
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
void aad::SetGradient(double g)
{
	var->grad = g;
}

void aad::Reset()
{
	double val = var->value;
	//delete var;
	var = new adjoint(val);
}

aad N(aad &f) // standard normal distribution
{
	aad res(0.5*(1 + erf(f.GetVar()->GetValue() / sqrt(2))));
	f.GetVar()->AppendToWeights(1.0 / sqrt(2.0 * M_PI)*exp(-0.5*pow(f.GetVar()->GetValue(), 2)));
	f.GetVar()->AppendToTape(res);
	return res;
}

aad sin(aad &f)
{
	aad res(sin(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(cos(f.GetVar()->GetValue()));
	f.GetVar()->AppendToTape(res);
	return res;
}

aad cos(aad &f)
{
	aad res(cos(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(-sin(f.GetVar()->GetValue()));
	f.GetVar()->AppendToTape(res);
	return res;
}

aad log(aad &f)
{
	aad res(log(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(1.0 / f.GetVar()->GetValue());
	f.GetVar()->AppendToTape(res);
	return res;
}

aad exp(aad &f)
{
	aad res(exp(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(exp(f.GetVar()->GetValue()));
	f.GetVar()->AppendToTape(res);
	return res;
}

aad sqrt(aad &f)
{
	aad res(sqrt(f.GetVar()->GetValue()));
	f.GetVar()->AppendToWeights(-0.5 / sqrt(f.GetVar()->GetValue()));
	f.GetVar()->AppendToTape(res);
	return res;
}







