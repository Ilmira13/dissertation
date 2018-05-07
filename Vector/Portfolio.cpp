#include <Portfolio.h>
vector<ForwardAD> VarVectorForwardAD(double &a, double &b, const int n)
{
	vector<ForwardAD> q(n);
	vector<double> d(n);
	for (int i = 0; i < n; i++)
	{
		q[i].value = a + (i + 1) * (b - a) / (n + 2);
		q[i].deriv = d;
		q[i].deriv[i] = 1; // derivative wrt a particular S[i]
	}	
	
	return  q;
}

vector<ReverseAD> VarVectorReverseAD(double &a, double &b, const int n)
{
	vector<ReverseAD> q(n);
	for (int i = 0; i < n; i++)
		q[i].GetVar()->SetValue(a + (i + 1) * (b - a) / (n + 2));

	return  q;
}

vector<double> VarVector(double &a, double &b,  const int n)
{
	vector<double> q(n);
	for (int i = 0; i < n; i++)
	{
		q[i] = a + (i + 1) * (b - a) / (n + 2);
	}
	return q;
}

vector<ForwardAD> portfolioAD(vector<ForwardAD> t, double &T, ForwardAD r, vector<ForwardAD> S, ForwardAD sigma, double &K)
{
	vector<ForwardAD> dd1 = d1(S, K, r, sigma, T, t);
	vector<ForwardAD> dd2 = d2(dd1, T, t, sigma);
	vector<ForwardAD> C = Call(dd1, dd2, S, K, r, sigma, T, t);
	return C;
}

vector<vector<double>> portfolioFD(vector<double> t, double &T, double &r, vector<double> S, double &sigma, double &K)
{
	vector<double> dd1 = d1(S, K, r, sigma, T, t);
	vector<double> dd2 = d2(dd1, T, t, sigma);
	vector<double> C = Call(dd1, dd2, S, K, r, sigma, T, t);

	vector<double> dd1S = d1(1.01*S, K, r, sigma, T, t);
	vector<double> dd2S = d2(dd1S, T, t, sigma);
	vector<double> CS = Call(dd1S, dd2S, 1.01*S, K, r, sigma, T, t);

	double sigma1 = 1.01*sigma;
	vector<double> dd1sigma = d1(S, K, r, sigma1, T, t);
	vector<double> dd2sigma = d2(dd1sigma, T, t, sigma1);
	vector<double> Csigma = Call(dd1sigma, dd2sigma, S, K, r, sigma1, T, t);

	vector<double> dd1t = d1(S, K, r, sigma, T, 1.01*t);
	vector<double> dd2t = d2(dd1t, T, 1.01*t, sigma);
	vector<double> Ct = Call(dd1t, dd2t, S, K, r, sigma, T, 1.01*t);

	double r1 = 1.01*r;
	vector<double> dd1r = d1(S, K, r1, sigma, T, t);
	vector<double> dd2r = d2(dd1r, T, t, sigma);
	vector<double> Cr = Call(dd1r, dd2r, S, K, r1, sigma, T, t);

	vector<vector<double>> f(S.size());
	for (int i = 0; i < S.size(); i++)
	{
		f[i].resize(5);
		f[i][0] = C[i]; // call
		f[i][1] = (CS[i] - C[i]) / (0.01*S[i]); // delta
		f[i][2] = (Csigma[i] - C[i]) / (0.01*sigma); // vega
		f[i][3] = (Ct[i] - C[i]) / (0.01*t[i]); // theta
		f[i][4] = (Cr[i] - C[i]) / (0.01*r); // rho
	}

	return f;
}


