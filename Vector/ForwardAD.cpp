#pragma hdrstop
#include <ForwardAD.h>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

vector<double> operator + (const vector<double> u)
{
	return u;
}

vector<double> operator - (const vector<double> u)
{
	vector<double> uu(size(u));
	for (int i = 0; i < size(u); i++)
	{
		uu[i] = -u[i];
	}
	return uu;
}

vector<double> operator + (const vector<double> u, const vector<double> h)
{
	int mini = min(u.size(), h.size());
	int maxi = max(u.size(), h.size());
	vector<double> F;
	if (maxi == u.size())
		F = u;
	else  F = h;
	vector<double> d(maxi);
	for (int i = 0; i < mini; i++)
	{
		d[i] = u[i] + h[i];
	}
	for (int i = mini; i < maxi; i++)
	{
		d[i] = F[i];
	}
	return d;
}

vector<double> operator - (const vector<double> u, const vector<double> h)
{
	int mini = min(u.size(), h.size());
	int maxi = max(u.size(), h.size());
	vector<double> F;
	if (maxi == u.size())
		F = u;
	else  F = h;
	vector<double> d(maxi);
	for (int i = 0; i < mini; i++)
	{
		d[i] = u[i] - h[i];
	}
	for (int i = mini; i < maxi; i++)
	{
		d[i] = F[i];
	}
	return d;
}

vector<double> operator * (const vector<double> u, const vector<double> h)
{
	int mini = min(u.size(), h.size());
	int maxi = max(u.size(), h.size());
	vector<double> F;
	if (maxi == u.size())
		F = u;
	else  F = h;
	vector<double> d(maxi);
	for (int i = 0; i < mini; i++)
	{
		d[i] = u[i] * h[i];
	}
	for (int i = mini; i < maxi; i++)
	{
		d[i] = F[i];
	}
	return d;
}

vector<double> operator * (const double &u, const vector<double> h)
{
	vector<double> k(h.size());
	for (int i = 0; i < h.size(); i++)
	{
		k[i] = u * h[i];
	}
	return k;
}

vector<double> operator * (const vector<double> h, const double &u)
{
	vector<double> k(h.size());
	for (int i = 0; i < h.size(); i++)
	{
		k[i] = u * h[i];
	}
	return k;
}

vector<double> operator / (const vector<double> u, const vector<double> h)
{
	int mini = min(u.size(), h.size());
	int maxi = max(u.size(), h.size());
	vector<double> F;
	if (maxi == u.size())
		F = u;
	else  F = h;
	vector<double> d(maxi);
	for (int i = 0; i < mini; i++)
	{
		d[i] = u[i] / h[i];
	}
	for (int i = mini; i < maxi; i++)
	{
		d[i] = F[i];
	}
	return d;
}

vector<double> operator / (const double &u, const vector<double> h)
{
	vector<double> k(h.size());
	for (int i = 0; i < h.size(); i++)
	{
		k[i] = u / h[i];
	}
	return k;
}

vector<double> operator / (const vector<double> h, const double &u)
{
	vector<double> k(h.size());
	for (int i = 0; i < h.size(); i++)
	{
		k[i] = h[i] / u;
	}
	return k;
}

ForwardAD::ForwardAD(double a, vector<double> b)
{
	value = a;
	deriv = b;
}

ForwardAD::ForwardAD()
{
	vector<double> f;
	value = 0;
	deriv = f;
}

ForwardAD::ForwardAD(double a)
{
	vector<double> f;
	value = a;
	deriv = f;
}

ForwardAD::ForwardAD(double val, int derivsCnt)
{
	vector<double> f(derivsCnt);
	value = val;
	deriv = f;
}


void ForwardAD::operator = (ForwardAD &c)
{
	this->value = c.value;
	this->deriv = c.deriv;
}

void ForwardAD::operator=(double c)
{
	this->value = c;
}



ForwardAD operator + (ForwardAD  &c)
{
	return c;
}

ForwardAD  operator - (const ForwardAD  &c)
{
	return ForwardAD(-c.value, -c.deriv);
}



ForwardAD  ForwardAD::operator * (const ForwardAD &c)
{
	return ForwardAD(this->value * c.value, this->deriv * c.value + this->value * c.deriv);
}
ForwardAD ForwardAD::operator * (const double &c)
{
	return ForwardAD(this->value * c, this->deriv * c);
}
ForwardAD operator * (const double c, const ForwardAD &f)
{
	return ForwardAD(c*f.value, c*f.deriv);
}

ForwardAD  ForwardAD::operator / (const ForwardAD &c)
{
	return ForwardAD(this->value / c.value, (this->deriv*c.value - this->value*c.deriv) / (c.value*c.value));
}
ForwardAD ForwardAD::operator / (const double &c)
{
	return ForwardAD(this->value / c, this->deriv / c);
}
bool ForwardAD::operator>(const ForwardAD & c) const
{
	return this->value > c.value;
}
bool ForwardAD::operator>=(const ForwardAD & c) const
{
	return this->value >= c.value;
}
bool ForwardAD::operator<(const ForwardAD & c) const
{
	return this->value < c.value;
}
bool ForwardAD::operator<=(const ForwardAD & c) const
{
	return this->value <= c.value;
}
bool ForwardAD::operator>(const double & c) const
{
	return this->value > c;
}
bool ForwardAD::operator>=(const double & c) const
{
	return this->value >= c;
}
bool ForwardAD::operator<(const double & c) const
{
	return this->value < c;
}
bool ForwardAD::operator<=(const double & c) const
{
	return this->value <= c;
}
ForwardAD operator / (const double c, const ForwardAD &f)
{
	return ForwardAD(c / f.value, -c * f.deriv / (f.value * f.value));
}

ForwardAD operator + (const ForwardAD &c, const ForwardAD &f)
{
	return ForwardAD(c.value + f.value, c.deriv + f.deriv);
}

ForwardAD operator + (const ForwardAD &f, const double c)
{
	return ForwardAD(c + f.value, f.deriv);
}

ForwardAD operator + (const double c, const ForwardAD &f)
{
	return ForwardAD(c + f.value, f.deriv);
}


ForwardAD operator - (const ForwardAD &c, const ForwardAD &f)
{
	return ForwardAD(c.value - f.value, c.deriv - f.deriv);
}
ForwardAD operator - (const ForwardAD &f, const double c)
{
	return ForwardAD(f.value - c, f.deriv);
}
ForwardAD operator - (const double c, const ForwardAD &f)
{
	return ForwardAD(c - f.value, -f.deriv);
}


ForwardAD N(ForwardAD &f) // standard normal distribution
{
	return ForwardAD(0.5*(1 + erf(f.value / sqrt(2))), dN(f.value)*f.deriv);
}

double N(const double &f) // standard normal distribution
{
	return 0.5*(1 + erf(f / sqrt(2)));
}

double dN(const double & f) // density
{
	return exp(-f * f / 2) / sqrt(2 * M_PI);
}

ForwardAD sin(ForwardAD &f)
{
	return ForwardAD(sin(f.value), cos(f.value)*f.deriv);
}

ForwardAD cos(ForwardAD &f)
{
	return ForwardAD(cos(f.value), -sin(f.value)*f.deriv);
}

ForwardAD log(ForwardAD &f)
{
	return ForwardAD(log(f.value), f.deriv / f.value);
}

ForwardAD exp(ForwardAD &f)
{
	return ForwardAD(exp(f.value), exp(f.value) * f.deriv);
}

ForwardAD sqrt(ForwardAD &f)
{
	return ForwardAD(sqrt(f.value), 1 / (2 * sqrt(f.value)) * f.deriv);
}
