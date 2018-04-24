#pragma once
#include <ForwardAD.h>
#include <ReverseAD.h>
#include <random> 

using namespace std;

enum MCtype
{
	MC1,
	MC2
};

void Reset(double &var);
void Reset(ForwardAD &var);
void Reset(ReverseAD &var);

double GetValue(double &var);
double GetValue(ForwardAD &var);
double GetValue(ReverseAD &var);

double GetGrad(double &var, int n = 0);
double GetGrad(ForwardAD &var, int n = 0);
double GetGrad(ReverseAD &var, int n = 0);

int GetVarsCount(double &var);
int GetVarsCount(ForwardAD &var);
int GetVarsCount(ReverseAD &var);

template <typename type1, typename type2, typename type3, typename type4, typename type5> class contract
{

private:
    default_random_engine generator;
    normal_distribution<double> distribution;

	type1* S;
	type2* sigma;
	type3* t;
	type4* r;
	type5 res;

	double price;
	double delta;
	double vega;
	double theta;
	double rho;
	vector<double> rnd;

	double K, T;
	type5 d1(); // we can do these functions private so that nobody could use them outside of this class
	type5 d2();
	type5 d1(type1 &St, type2 &sigmat, type3 &tt, type4 &rt);
	type5 d2(type1 &St, type2 &sigmat, type3 &tt, type4 &rt);

	type5 MCPath1(type1 &St, type2 &sigmat, type3 &rt, type4 & tt, double rnd, int n);
	type5 MCPath2(type1 &St, type2 &sigmat, type3 &rt, type4 & tt, double rnd, int n);

public:
	contract(type1*, type2*, type3*, type4*, double, double);
	contract();
	void CalculateGradient();
	void CalculatePrice();
	void CalculatePriceMC(int, int, MCtype, bool useTheSameRndSequence=false);
	double GetPrice();
	double GetDelta();
	double GetVega();
	double GetRho();
	double GetTheta();

	type5 Delta();
	type5 Vega();
	type5 Theta();
	type5 Rho();

};

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline contract<type1, type2, type3, type4, type5>::contract(type1* S1, type2* sigma1, type3* t1, type4* r1, double T1, double K1):
	price(0.0)
{
	S = S1; sigma = sigma1; t = t1;	r = r1;	T = T1;	K = K1;
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline contract<type1, type2, type3, type4, type5>::contract() : distribution(0.0, 1.0), price(0.0)
{
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline void contract<type1, type2, type3, type4, type5>::CalculateGradient()
{
	//if (typeid(type5) == typeid(ForwardAD))
	//{
		delta = 0.0;
		vega = 0.0;
		theta = 0.0;
		rho = 0.0;
		for (int i = 0; i < GetVarsCount(res); ++i)
		{
			delta += GetGrad(res, i) * GetGrad(*S, i);
			vega += GetGrad(res, i) * GetGrad(*sigma, i);
			theta += GetGrad(res, i) * GetGrad(*t, i);
			rho += GetGrad(res, i) * GetGrad(*r, i);
		}
	/*}
	else if (typeid(type5) == typeid(ReverseAD))
	{
		delta = GetGrad(*S);
		vega = GetGrad(*sigma);
		theta = GetGrad(*t);
		rho = GetGrad(*r);
	}
	else if (typeid(type5) == typeid(double))
	{
		delta = 0.0;
		vega = 0.0;
		theta = 0.0;
		rho = 0.0;
	}
	else
	{
		throw exception::exception("Unknown variable type");
	}*/
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline void contract<type1, type2, type3, type4, type5>::CalculatePrice()
{
	res = *S * N(d1()) - K * exp(-*r * (T - *t))*N(d2());
	price = GetValue(res);
	CalculateGradient();
}


template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::d1()
{
	return (log(*S / K) + (*r + *sigma*(*sigma) / 2.0)*(T - *t)) / (*sigma*sqrt(T - *t));
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::d1(type1 &St, type2 &sigmat, type3 &tt, type4 &rt)
{
	return (log(St / K) + (rt + sigmat*sigmat / 2.0)*(T - tt)) / (sigmat*sqrt(T - tt));
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::d2()
{
	return (d1() - *sigma*sqrt(T - *t));
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::d2(type1 &St, type2 &sigmat, type3 &tt, type4 &rt)
{
	return (d1(St, sigmat, tt, rt) - sigmat*sqrt(T - tt));
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
inline type5 contract<type1, type2, type3, type4, type5>::MCPath1(type1 &St, type2 &sigmat, type3 &rt, type4 & tt, double rnd, int n)
{
	return St * (1 + rt * (T - tt) / n + sigmat * sqrt((T - tt) / n)*rnd);
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline type5 contract<type1, type2, type3, type4, type5>::MCPath2(type1 &St, type2 &sigmat, type3 &rt, type4 & tt, double rnd, int n)
{
	return St * exp((rt - sigmat*sigmat / 2)*(T - tt) / n + sigmat*sqrt((T - tt) / n) * rnd);
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline void contract<type1, type2, type3, type4, type5>::CalculatePriceMC(int n, int m, MCtype type, bool useTheSameRndSequence)
{
	// n - time
	// m - trajectory
	/*
	price = 0.0; // potentially may cause problems in derivatives calculation with ForwardAD
	delta = 0.0;
	vega = 0.0;
	rho = 0.0;
	theta = 0.0;
	vector<type5> result(m + 1, 0.0);
	vector<type1> St(n, *S);
	type2 sigmat(*sigma);
	type3 tt(*t);
	type4 rt(*r);

	if (!useTheSameRndSequence || size(rnd) == 0)
	{
		rnd.resize(m*(n - 1));
		for (int i = 0; i < m; i++)
			for (int j = 0; j < (n - 1); j++)
				rnd[i*(n - 1) + j] = distribution(generator);
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 1; j < n; j++)
		{
			switch (type)
			{
			case MCtype::MC1:
				St[j] = MCPath1(St[j - 1], sigmat, rt, tt, rnd[i*(n - 1) + j - 1], n);
				break;
			case MCtype::MC2:
				St[j] = MCPath2(St[j - 1], sigmat, rt, tt, rnd[i*(n - 1) + j - 1], n);
				break;
			default:
				throw exception("unknown MC type: " + type);
				break;
			}
		}
		
		result[i + 1] = result[i] + exp(-rt*(T - tt))*(St[n - 1] - K)*((St[n - 1] > K) ? 1.0 : 0.0);
		//if (St[n - 1] > K)
			//res[i + 1] = res[i] + exp(-rt * (T - tt))*(St[n - 1] - K);
	}

	res = result[m] / m;
	price = GetValue(res);

	delta = GetGrad(St[0]) / m;
	vega = GetGrad(sigmat) / m;
	rho = GetGrad(rt) / m;
	theta = GetGrad(tt) / m;
	*/
	vector<type5> result(m + 1, 0.0);
	vector<type1> St(n);
	type2 sigmat;
	type3 tt;
	type4 rt;

	St[0] = *S;
	sigmat = *sigma;
	tt = *t;
	rt = *r;

	if (!useTheSameRndSequence || size(rnd) == 0)
	{
		rnd.resize(m*(n - 1));
		for (int i = 0; i < m; i++)
			for (int j = 0; j < (n - 1); j++)
				rnd[i*(n - 1) + j] = distribution(generator);
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 1; j < n; j++)
		{
			switch (type)
			{
			case MCtype::MC1:
				St[j] = MCPath1(St[j - 1], sigmat, rt, tt, rnd[i*(n - 1) + j - 1], n);
				break;
			case MCtype::MC2:
				St[j] = MCPath2(St[j - 1], sigmat, rt, tt, rnd[i*(n - 1) + j - 1], n);
				break;
			default:
				throw exception("unknown MC type: " + type);
				break;
			}
		}

		result[i + 1] = result[i] + exp(-rt * (T - tt))*(St[n - 1] - K)*((St[n - 1] > K) ? 1.0 : 0.0);
	}

	res = result[m] / m;
	price = GetValue(res);
	CalculateGradient();
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline double contract<type1, type2, type3, type4, type5>::GetPrice()
{
	return price;
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline double contract<type1, type2, type3, type4, type5>::GetDelta()
{
	return delta;
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline double contract<type1, type2, type3, type4, type5>::GetVega()
{
	return vega;
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline double contract<type1, type2, type3, type4, type5>::GetRho()
{
	return rho;
}

template<typename type1, typename type2, typename type3, typename type4, typename type5>
inline double contract<type1, type2, type3, type4, type5>::GetTheta()
{
	return theta;
}



vector<ForwardAD> d1(vector<ForwardAD> S, double &K, ForwardAD r, ForwardAD sigma, double &T, vector<ForwardAD> t);
vector<ForwardAD> d2(vector<ForwardAD> d1, double &T, vector<ForwardAD> t, ForwardAD sigma);
vector<ForwardAD> Call(vector<ForwardAD> d1, vector<ForwardAD> d2, vector<ForwardAD> S, double &K, ForwardAD r, ForwardAD sigma, double &T, vector<ForwardAD> t);


vector<double> d1(vector<double> S, double &K, double &r, double &sigma, double &T, vector<double> t);
vector<double> d2(vector<double> d1, double &T, vector<double> t, double &sigma);
vector<double> Call(vector<double> d1, vector<double> d2, vector<double> S, double &K, double &r, double &sigma, double &T, vector<double> t);

