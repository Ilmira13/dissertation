#pragma once
#include <ForwardAD.h>
#include <ReverseAD.h>
#include <BlackScholes.h>

vector<ForwardAD> VarVectorForwardAD(double &a, double &b, const int n);
vector<ReverseAD> VarVectorReverseAD(double &a, double &b, const int n);
vector<double> VarVector(double &a, double &b, const int n);
vector<ForwardAD> portfolioAD(vector<ForwardAD> t, double &T, ForwardAD r, vector<ForwardAD> S, ForwardAD sigma, double &K);
vector<vector<double>> portfolioFD(vector<double> t, double &T, double &r, vector<double> S, double &sigma, double &K);

using namespace std;
template<typename contracttype, typename resulttype> class portfolio
{
private:
	double price;
	vector<double> deltas;
	vector<double> vegas;
	vector<double> rhos;
	vector<double> thetas;
	int j; // time
	int m; // trajectory
	string PricingMethod;
	vector<contracttype> contracts;
	void CalculatePPrice();
	void CalculatePPriceMC1(int, int, bool useTheSameRndSequence = false);
	void CalculatePPriceMC2(int, int, bool useTheSameRndSequence = false);

public:
	portfolio();
	portfolio(vector<contracttype> &, string);
	portfolio(vector<contracttype> &, string, int, int);
	void Price(bool useTheSameRndSequence = false);
	vector<double> FiniteDiff (vector<double> &, int);
	vector<resulttype> Delta();
	vector<resulttype> Vega();
	vector<resulttype> Theta();
	vector<resulttype> Rho();

	double GetPrice();
	vector<double> GetDeltas();
	vector<double> GetVegas();
	vector<double> GetRhos();
	vector<double> GetThetas();
};

template<typename contracttype, typename resulttype> inline portfolio<contracttype, resulttype>::portfolio() 
{
	
}

template<typename contracttype, typename resulttype> inline portfolio<contracttype, resulttype>::portfolio(vector<contracttype> &iContracts, string str)
{
	contracts = iContracts;
	PricingMethod = str;
	deltas.resize(size(contracts));
	vegas.resize(size(contracts));
	rhos.resize(size(contracts));
	thetas.resize(size(contracts));
}

template<typename contracttype, typename resulttype>
inline portfolio<contracttype, resulttype>::portfolio(vector<contracttype>&iContracts, string str, int jj, int mm)
{
	contracts = iContracts;
	PricingMethod = str;
	deltas.resize(size(contracts));
	vegas.resize(size(contracts));
	rhos.resize(size(contracts));
	thetas.resize(size(contracts));
	j = jj;
	m = mm;
}

template<typename contracttype, typename resulttype>
inline void portfolio<contracttype, resulttype>::Price(bool useTheSameRndSequence)
{
	if (PricingMethod == "MC1")
	{
		this->CalculatePPriceMC1(j, m, useTheSameRndSequence);
	}
	else if (PricingMethod == "MC2")
	{
		this->CalculatePPriceMC2(j, m, useTheSameRndSequence);
	}
	else 
		this->CalculatePPrice();
}


template<typename contracttype, typename resulttype>
inline void portfolio<contracttype, resulttype>::CalculatePPrice()
{
	price = 0.0;
	for (int i = 0; i < size(contracts); ++i)
	{
		contracts[i].CalculatePrice();
		price = price + contracts[i].GetPrice();
		deltas[i] = contracts[i].GetDelta();
		vegas[i] = contracts[i].GetVega();
		rhos[i] = contracts[i].GetRho();
		thetas[i] = contracts[i].GetTheta();
	}
}

template<typename contracttype, typename resulttype>
inline void portfolio<contracttype, resulttype>::CalculatePPriceMC1(int j, int m, bool useTheSameRndSequence)
{
	price = 0.0;
	for (int i = 0; i < size(contracts); ++i)
	{
		contracts[i].CalculatePriceMC(j, m, MCtype::MC1, useTheSameRndSequence);
		price = price + contracts[i].GetPrice();
		deltas[i] = contracts[i].GetDelta();
		vegas[i] = contracts[i].GetVega();
		rhos[i] = contracts[i].GetRho();
		thetas[i] = contracts[i].GetTheta();
	}
}

template<typename contracttype, typename resulttype>
inline void portfolio<contracttype, resulttype>::CalculatePPriceMC2(int j, int m, bool useTheSameRndSequence)
{
	price = 0.0;
	for (int i = 0; i < size(contracts); ++i)
	{
		contracts[i].CalculatePriceMC(j, m, MCtype::MC2, useTheSameRndSequence);
		price = price + contracts[i].GetPrice();
		deltas[i] = contracts[i].GetDelta();
		vegas[i] = contracts[i].GetVega();
		rhos[i] = contracts[i].GetRho();
		thetas[i] = contracts[i].GetTheta();
	}
}


template<typename contracttype, typename resulttype>
inline vector<double> portfolio<contracttype, resulttype>::FiniteDiff(vector<double> &var, int bmp)
{
	this->Price(true);
	double tmp = 0;
	double bump = double(bmp) / 100;
	double priceBase = 0.0;

	vector<double> derivs(size(var), 0.0);
	for (int i = 0; i < size(var); ++i)
	{
		tmp = var[i];
		priceBase = price;
		var[i] = var[i] * (1 + bump);
		this->Price(true);
		derivs[i] = (price - priceBase) / (var[i] - tmp); // for MC we should use the same seed!
		var[i] = tmp;
		price = priceBase;
	}
	return derivs;
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

template<typename contracttype, typename resulttype>
inline double portfolio<contracttype, resulttype>::GetPrice()
{
	return price;
}

template<typename contracttype, typename resulttype>
inline vector<double> portfolio<contracttype, resulttype>::GetDeltas()
{
	return deltas;
}

template<typename contracttype, typename resulttype>
inline vector<double> portfolio<contracttype, resulttype>::GetVegas()
{
	return vegas;
}

template<typename contracttype, typename resulttype>
inline vector<double> portfolio<contracttype, resulttype>::GetRhos()
{
	return rhos;
}

template<typename contracttype, typename resulttype>
inline vector<double> portfolio<contracttype, resulttype>::GetThetas()
{
	return thetas;
}
