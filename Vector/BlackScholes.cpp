#include <BlackScholes.h>

void Reset(double &var)
{
}

void Reset(ForwardAD &var)
{
}

void Reset(ReverseAD &var)
{
	var.Reset();
}

double GetValue(double &var)
{
	return var;
}

double GetValue(ForwardAD &var)
{
	return var.value;
}

double GetValue(ReverseAD &var)
{
	return var.GetVar()->GetValue();
}

double GetGrad(double &var, int n)
{
	return 0.0;
}

double GetGrad(ForwardAD &var, int n)
{
	return var.deriv[n];
}

double GetGrad(ReverseAD &var, int n)
{
	return var.GetGradient();
}

int GetVarsCount(double &var)
{
	return 1.0;
}

int GetVarsCount(ForwardAD &var)
{
	return size(var.deriv);
}

int GetVarsCount(ReverseAD &var)
{
	return 1.0;
}

vector<ForwardAD> d1(vector<ForwardAD> S, double &K, ForwardAD r, ForwardAD sigma, double &T, vector<ForwardAD> t)
{
	vector<ForwardAD> f(S.size());
	for (int i = 0; i < S.size(); i++)
	{
		f[i] = (log(S[i] / K) + (r + sigma*sigma / 2)*(T - t[i])) / (sigma*sqrt(T - t[i]));
	}
	return f;
}

vector<ForwardAD> d2(vector<ForwardAD> d1, double &T, vector<ForwardAD>t, ForwardAD sigma)
{
	vector<ForwardAD> f(d1.size());
	for (int i = 0; i < d1.size(); i++)
	{
		f[i] = d1[i] - sigma*sqrt(T - t[i]);
	}
	return f;
}

vector<ForwardAD> Call(vector<ForwardAD> d1, vector<ForwardAD> d2, vector<ForwardAD> S, double &K, ForwardAD r, ForwardAD sigma, double &T, vector<ForwardAD> t)
{
	vector<ForwardAD> f(d1.size());
	for (int i = 0; i < d1.size(); i++)
	{
		f[i] = S[i]*N(d1[i]) - K*exp(-r*(T - t[i]))*N(d2[i]);
	}
	return f;
}

vector<double> d1(vector<double> S, double &K, double &r, double &sigma, double &T, vector<double> t)
{
	vector<double> f(S.size());
	for (int i = 0; i < S.size(); i++)
	{
		f[i] = (log(S[i] / K) + (r + sigma*sigma / 2)*(T - t[i])) / (sigma*sqrt(T - t[i]));
	}
	return f;
}

vector<double> d2(vector<double> d1, double &T, vector<double> t, double &sigma)
{
	vector<double> f(d1.size());
	for (int i = 0; i < d1.size(); i++)
	{
		f[i] = d1[i] - sigma*sqrt(T - t[i]);
	}
	return f;
}

vector<double> Call(vector<double> d1, vector<double> d2, vector<double> S, double &K, double &r, double &sigma, double &T, vector<double> t)
{
	vector<double> f(d1.size());
	for (int i = 0; i < d1.size(); i++)
	{
		f[i] = S[i] * N(d1[i]) - K*exp(-r*(T - t[i]))*N(d2[i]);
	}
	return f;
}