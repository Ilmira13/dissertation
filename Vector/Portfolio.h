#pragma once
#include <autodiff.h>
#include <AAD.h>
#include <BlackScholes.h>

vector<autodiff> VarVectorAD(double &a, double &b, const int n);
vector<aad> VarVectorAAD(double &a, double &b, const int n);
vector<double> VarVector(double &a, double &b, const int n);
vector<autodiff> portfolioAD(vector<autodiff> t, double &T, autodiff r, vector<autodiff> S, autodiff sigma, double &K);
vector<vector<double>> portfolioFD(vector<double> t, double &T, double &r, vector<double> S, double &sigma, double &K);

using namespace std;
template<typename contracttype, typename resulttype> class portfolio
{
private:
	int j; // time
	int m; // trajectory
	string PricingMethod;
	vector<contracttype> contracts;
	resulttype CalculatePPrice();
	resulttype CalculatePPriceMC1(int, int);
	resulttype CalculatePPriceMC2(int, int);

public:
	portfolio();
	portfolio(vector<contracttype> &, string);
	portfolio(vector<contracttype> &, string, int, int);
	resulttype Price();
	vector<resulttype> FiniteDiff (double, vector<double> &, int);
	vector<resulttype> Delta();
	vector<resulttype> Vega();
	vector<resulttype> Theta();
	vector<resulttype> Rho();
};

template<typename contracttype, typename resulttype> inline portfolio<contracttype, resulttype>::portfolio() 
{
	
}

template<typename contracttype, typename resulttype> inline portfolio<contracttype, resulttype>::portfolio(vector<contracttype> &iContracts, string str)
{
	contracts = iContracts;
	PricingMethod = str;
}

template<typename contracttype, typename resulttype>
inline portfolio<contracttype, resulttype>::portfolio(vector<contracttype>&iContracts, string str, int jj, int mm)
{
	contracts = iContracts;
	PricingMethod = str;
	j = jj;
	m = mm;
}

template<typename contracttype, typename resulttype>
inline resulttype portfolio<contracttype, resulttype>::Price()
{
	if (PricingMethod == "MC1")
	{
		return this->CalculatePPriceMC1(j, m);
	}
	else if (PricingMethod == "MC2")
	{
		return this->CalculatePPriceMC2(j, m);
	}
	else 
		return this->CalculatePPrice();
}


template<typename contracttype, typename resulttype>
inline resulttype portfolio<contracttype, resulttype>::CalculatePPrice()
{
	resulttype res = (resulttype)contracts[0].CalculatePrice();
	for (int i = 1; i < size(contracts); ++i)
	{
		res = res + (resulttype)contracts[i].CalculatePrice();
	}
	return res;
}

template<typename contracttype, typename resulttype>
inline resulttype portfolio<contracttype, resulttype>::CalculatePPriceMC1(int j, int m)
{
	resulttype res = (resulttype)contracts[0].CalculatePriceMC1(j, m);
	for (int i = 1; i < size(contracts); ++i)
	{
		res = res + (resulttype)contracts[i].CalculatePriceMC1(j, m);
	}
	return res;
}

template<typename contracttype, typename resulttype>
inline resulttype portfolio<contracttype, resulttype>::CalculatePPriceMC2(int j, int m)
{
	resulttype res = (resulttype)contracts[0].CalculatePriceMC2(j, m);
	for (int i = 1; i < size(contracts); ++i)
	{
		res = res + (resulttype)contracts[i].CalculatePriceMC2(j, m);
	}
	return res;
}


template<typename contracttype, typename resulttype>
inline vector<resulttype> portfolio<contracttype, resulttype>::FiniteDiff(double priceDBL, vector<double> &var, int bmp)
{
	vector<double> derivsFD(size(contracts));
	double tmp = 0;
	double bump = double(bmp)/100;

		for (int i = 0; i < size(contracts); ++i)
		{
			tmp = var[i];
			var[i] = var[i] * (1 + bump);
			derivsFD[i] = (this->Price() - priceDBL) / (var[i] - tmp);
			var[i] = tmp;
		}
	return derivsFD;
}

template<typename contracttype, typename resulttype>
inline vector<resulttype> portfolio<contracttype, resulttype>::Delta()
{
	vector<resulttype> delta(size(contracts));
	for (int i = 0; i < size(contracts); i++)
	{
		delta[i] = contracts[i].Delta();
	}
	return delta;
}


template<typename contracttype, typename resulttype>
inline vector<resulttype> portfolio<contracttype, resulttype>::Vega()
{
	vector<resulttype> delta(size(contracts));
	for (int i = 0; i < size(contracts); i++)
	{
		delta[i] = contracts[i].Vega();
	}
	return delta;
}

template<typename contracttype, typename resulttype>
inline vector<resulttype> portfolio<contracttype, resulttype>::Theta()
{
	vector<resulttype> delta(size(contracts));
	for (int i = 0; i < size(contracts); i++)
	{
		delta[i] = contracts[i].Theta();
	}
	return delta;
}

template<typename contracttype, typename resulttype>
inline vector<resulttype> portfolio<contracttype, resulttype>::Rho()
{
	vector<resulttype> delta(size(contracts));
	for (int i = 0; i < size(contracts); i++)
	{
		delta[i] = contracts[i].Rho();
	}
	return delta;
}
