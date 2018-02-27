#include <BlackScholes.h>





vector<autodiff> d1(vector<autodiff> S, double &K, autodiff r, autodiff sigma, double &T, vector<autodiff> t)
{
	vector<autodiff> f(S.size());
	for (int i = 0; i < S.size(); i++)
	{
		f[i] = (log(S[i] / K) + (r + sigma*sigma / 2)*(T - t[i])) / (sigma*sqrt(T - t[i]));
	}
	return f;
}

vector<autodiff> d2(vector<autodiff> d1, double &T, vector<autodiff>t, autodiff sigma)
{
	vector<autodiff> f(d1.size());
	for (int i = 0; i < d1.size(); i++)
	{
		f[i] = d1[i] - sigma*sqrt(T - t[i]);
	}
	return f;
}

vector<autodiff> Call(vector<autodiff> d1, vector<autodiff> d2, vector<autodiff> S, double &K, autodiff r, autodiff sigma, double &T, vector<autodiff> t)
{
	vector<autodiff> f(d1.size());
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