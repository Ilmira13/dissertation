#pragma once
#include <BlackScholes.h>
vector<double> delta (vector<double> d1)
{
	vector<double> f(d1.size());
	for (int i = 0; i < d1.size(); i++)
	{
		f[i] = N(d1[i]); // delta (S) 	
	}
	return f;
}

vector<double> vega (vector<double> S, vector<double> d1, double &T, vector<double> t)
{
	vector<double> f(S.size());
	for (int i = 0; i < S.size(); i++)
	{
		f[i] = S[i] * exp(-d1[i] * d1[i] / 2) / sqrt(2 * M_PI) * sqrt(T - t[i]); // vega (sigma);
	}
	return f;
}

vector<double> theta (vector<double> S, vector<double> d1, double &sigma, double &T, vector<double> t, double &r, double &K, vector<double> d2)
{
	vector<double> f(S.size());
	for (int i = 0; i < S.size(); i++)
	{
		f[i] = -S[i] * exp(-0.5*d1[i] * d1[i]) * sigma / (2 * sqrt(2 * M_PI*(T - t[i]))) - r*K*exp(-r*(T - t[i]))*N(d2[i]); // theta (t)
	}
	return f;
}

vector<double> rho (double &K, double &T, double &r, vector<double> t, vector<double> d2)
{
	vector<double> f(d2.size());
	for (int i = 0; i < d2.size(); i++)
	{
		f[i] = K*(T - t[i]) * exp(-r*(T - t[i]))*N(d2[i]); // rho (r)
	}
	return f;
}