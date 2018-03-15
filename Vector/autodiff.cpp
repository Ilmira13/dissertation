#pragma hdrstop
#include <autodiff.h>
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
	vector<double> uu(u.size());
	for (int i = 0; i < u.size(); i++)
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

autodiff::autodiff(double a, vector<double> b)
{
	value = a;
	deriv = b;
}

autodiff::autodiff()
{
	vector<double> f;
	value = 0;
	deriv = f;
}

autodiff::autodiff(double a)
{
	vector<double> f;
	value = a;
	deriv = f;
}

autodiff::autodiff(double val, int derivsCnt)
{
	vector<double> f(derivsCnt);
	value = val;
	deriv = f;
}


void autodiff::operator = (autodiff &c)
{
	this->value = c.value;
	this->deriv = c.deriv;
}

void autodiff::operator=(double c)
{
	this->value = c;
}



autodiff operator + (autodiff  &c)
{
	return c;
}

autodiff  operator - (const autodiff  &c)
{
	return autodiff(-c.value, -c.deriv);
}



autodiff  autodiff::operator * (const autodiff &c)
{
	return autodiff(this->value * c.value, this->deriv * c.value + this->value * c.deriv);
}
autodiff autodiff::operator * (const double &c)
{
	return autodiff(this->value * c, this->deriv * c);
}
autodiff operator * (const double c, const autodiff &f)
{
	return autodiff(c*f.value, c*f.deriv);
}

autodiff  autodiff::operator / (const autodiff &c)
{
	return autodiff(this->value / c.value, (this->deriv*c.value - this->value*c.deriv) / (c.value*c.value));
}
autodiff autodiff::operator / (const double &c)
{
	return autodiff(this->value / c, this->deriv / c);
}
bool autodiff::operator>(const autodiff & c) const
{
	return this->value > c.value;
}
bool autodiff::operator>=(const autodiff & c) const
{
	return this->value >= c.value;
}
bool autodiff::operator<(const autodiff & c) const
{
	return this->value < c.value;
}
bool autodiff::operator<=(const autodiff & c) const
{
	return this->value <= c.value;
}
bool autodiff::operator>(const double & c) const
{
	return this->value > c;
}
bool autodiff::operator>=(const double & c) const
{
	return this->value >= c;
}
bool autodiff::operator<(const double & c) const
{
	return this->value < c;
}
bool autodiff::operator<=(const double & c) const
{
	return this->value <= c;
}
autodiff operator / (const double c, const autodiff &f)
{
	return autodiff(c / f.value, -c * f.deriv / (f.value * f.value));
}

autodiff operator + (const autodiff &c, const autodiff &f)
{
	return autodiff(c.value + f.value, c.deriv + f.deriv);
}

autodiff operator + (const autodiff &f, const double c)
{
	return autodiff(c + f.value, f.deriv);
}

autodiff operator + (const double c, const autodiff &f)
{
	return autodiff(c + f.value, f.deriv);
}


autodiff operator - (const autodiff &c, const autodiff &f)
{
	return autodiff(c.value - f.value, c.deriv - f.deriv);
}
autodiff operator - (const autodiff &f, const double c)
{
	return autodiff(f.value - c, f.deriv);
}
autodiff operator - (const double c, const autodiff &f)
{
	return autodiff(c - f.value, -f.deriv);
}


autodiff N(autodiff &f) // standard normal distribution
{
	return autodiff(0.5*(1 + erf(f.value / sqrt(2))), dN(f.value)*f.deriv);
}

double N(const double &f) // standard normal distribution
{
	return 0.5*(1 + erf(f / sqrt(2)));
}

double dN(const double & f) // density
{
	return exp(-f * f / 2) / sqrt(2 * M_PI);
}

autodiff sin(autodiff &f)
{
	return autodiff(sin(f.value), cos(f.value)*f.deriv);
}

autodiff cos(autodiff &f)
{
	return autodiff(cos(f.value), -sin(f.value)*f.deriv);
}

autodiff log(autodiff &f)
{
	return autodiff(log(f.value), f.deriv / f.value);
}

autodiff exp(autodiff &f)
{
	return autodiff(exp(f.value), exp(f.value) * f.deriv);
}

autodiff sqrt(autodiff &f)
{
	return autodiff(sqrt(f.value), 1 / (2 * sqrt(f.value)) * f.deriv);
}















