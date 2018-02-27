#pragma once
#include <autodiff.h>
#include <random> 


using namespace std;

template <typename type1, typename type2, typename type3, typename type4, typename type5> class contract
{

private:
    default_random_engine generator;
    normal_distribution<double> distribution;
	type1* S;
	type2* sigma;
	type3* t;
	type4* r;
	double K, T;
	type5 d1(); // we can do these functions private so that nobody could use them outside of this class
	type5 d2();

public:
	contract(type1*, type2*, type3*, type4*, double, double);
	contract();
	type5 CalculatePrice();
	type5 CalculatePriceMC1(int, int);
	type5 CalculatePriceMC2(int, int);

	type5 Delta();
	type5 Vega();
	type5 Theta();
	type5 Rho();

};

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline contract<type1, type2, type3, type4, type5>::contract(type1* S1, type2* sigma1, type3* t1, type4* r1, double T1, double K1)
{
	S = S1; sigma = sigma1; t = t1;	r = r1;	T = T1;	K = K1;
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline contract<type1, type2, type3, type4, type5>::contract() : distribution(0.0, 1.0)
{
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::CalculatePrice()
{
	return *S * N(d1()) - K*exp(-*r*(T - *t))*N(d2());
}


template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::d1()
{
	return (log(*S / K) + (*r + *sigma*(*sigma) / 2.0)*(T - *t)) / (*sigma*sqrt(T - *t));
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::d2()
{
	return (d1() - *sigma*sqrt(T - *t));
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::Delta()
{
	return N(d1());
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::Vega()
{
	return dN(d1()) **S*sqrt(T - *t);
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::Theta()
{
	return -*S*dN(d1())**sigma/2/sqrt(T - *t) - *r*K*N(d2())*exp(-*r*(T - *t));
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::Rho()
{
	return K*(T - *t)*N(d2())*exp(-*r*(T - *t));
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::CalculatePriceMC1(int n, int m)
{
	// n - time
	// m - trajectory

	type5 res(0.0); // potentially may cause problems in derivatives calculation with autodiff
	for (int i = 0; i < m; i++)
	{
		type1 St = *S;
		for (int j = 1; j < n; j++)
		{
			St = St * (1 + *r*(T / n) + *sigma*sqrt((T / n))*distribution(generator));
		}
		if (St > K)
			res = res + (St - K);
	}
	res = res / m;
	res = exp(-*r*(T - *t))*res;
	return res;
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::CalculatePriceMC2(int n, int m)
{
	// n - time
	// m - trajectory

	type5 res(0.0); // potentially may cause problems in derivatives calculation with autodiff
	for (int i = 0; i < m; i++)
	{
		type1 St = *S;
		for (int j = 1; j < n; j++)
		{
			St = St * exp((*r - *sigma**sigma / 2)*(T / n) + *sigma*sqrt(T / n) * distribution(generator));
		}
		if (St > K)
			res = res + (St - K);
	}
	res = res / m;
	res = exp(-*r*(T - *t))*res;
	return res;
}



vector<autodiff> d1(vector<autodiff> S, double &K, autodiff r, autodiff sigma, double &T, vector<autodiff> t);
vector<autodiff> d2(vector<autodiff> d1, double &T, vector<autodiff> t, autodiff sigma);
vector<autodiff> Call(vector<autodiff> d1, vector<autodiff> d2, vector<autodiff> S, double &K, autodiff r, autodiff sigma, double &T, vector<autodiff> t);


vector<double> d1(vector<double> S, double &K, double &r, double &sigma, double &T, vector<double> t);
vector<double> d2(vector<double> d1, double &T, vector<double> t, double &sigma);
vector<double> Call(vector<double> d1, vector<double> d2, vector<double> S, double &K, double &r, double &sigma, double &T, vector<double> t);

