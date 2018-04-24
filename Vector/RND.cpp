#pragma once
#include <RND.h>

void SpotTest(const int it, const int n)
{
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, n);
	vector<ForwardAD> SAD = VarVectorAD(aS, bS, n); // will not work if we have other AD variables except S
	double K = 300;
	double Sigma = 0.25;
	double r = 0.03;
	double T = 1;
	double t = 1. / 365;

	clock_t time = clock();

	// compose 2 portfolios: one for double arguments and one for ForwardAD arguments
	vector<contract<double, double, double, double, double> > contractsDBL;
	vector<contract<ForwardAD, double, double, double, ForwardAD> > contractsAD;
	for (int i = 0; i < n; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<ForwardAD, double, double, double, ForwardAD> cAD(&SAD[i], &Sigma, &t, &r, T, K);
		contractsAD.push_back(cAD);
	}
	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	portfolio<contract<ForwardAD, double, double, double, ForwardAD>, ForwardAD> P_AD(contractsAD, "0");
	cout << "AD function calculation time: ";
	time = clock();
	for (int i = 0; i < it; ++i)
		P_AD.Price();
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "DBL function calculation time: ";
	time = clock();
	for (int i = 0; i < it; ++i)
		P_DBL.Price();
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

	// finite-difference derivatives
	vector<double> derivsFD(n);
	cout << "FD derivatives calculation time: ";
	time = clock();
	for (int i = 0; i < it; ++i)
		derivsFD = P_DBL.FiniteDiff(S, 1);
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

	cout << "DBL function value = " << P_DBL.GetPrice() << endl;
	cout << "AD function value = " << P_AD.GetPrice() << endl;


	// printing out derivative values: AD versus analytical vs FD
	double totalDelta = 0;
	cout << "Analytical deltas:" << endl;
	for (int i = 0; i < size(contractsDBL); ++i)
	{
		totalDelta += contractsDBL[i].Delta();
		cout << contractsDBL[i].Delta() << "\t";
	}
	cout << endl << "Portfolio analytical total delta = " << totalDelta << endl;

	totalDelta = 0;
	cout << "AD deltas:" << endl;
	vector<double> deltas = P_AD.GetDeltas();
	for (int i = 0; i < n; ++i) // in the case of this test we have n independent variables
	{
		totalDelta += deltas[i];
		cout << deltas[i] << "\t";
	}
	cout << endl << "Portfolio AD total delta = " << totalDelta << endl;

	totalDelta = 0;
	cout << "FD deltas:" << endl;
	for (int i = 0; i < n; ++i) // in the case of this test we have n independent variables
	{
		totalDelta += derivsFD[i];
		cout << derivsFD[i] << "\t";
	}
	cout << endl << "Portfolio AD total delta = " << totalDelta << endl;

	cout << deltas[0] << endl;
	cout << deltas[2] << endl;
}

void AllVarTest(const int it, const int a)
{
	
	double aS = 100; double bS = 500;
	double aSigma = 0.15; double bSigma = 0.3;
	double aR = 0.01; double bR = 0.04;
	double aT = 1. / 365;
	double K = 300;
	double T = 1;
	vector<double> S = VarVector(aS, bS, a); vector<ForwardAD>  SAD = VarVectorAD(aS, bS, a);
	vector<double> sigma = VarVector(aSigma, bSigma, a); vector<ForwardAD> sigmReverseAD = VarVectorAD(aSigma, bSigma, a);
	vector<double> 	r = VarVector(aR, bR, a); vector<ForwardAD>  rAD = VarVectorAD(aR, bR, a);
	vector<double> t = VarVector(aT, aT, a); vector<ForwardAD> tAD = VarVectorAD(aT, aT, a);
	vector<vector<double>> derivs(4 * a);
	for (int i = 0; i < size(derivs); i++)
	{
		derivs[i].resize(4 * a);
		for (int j = 0; j < size(derivs); j++)
		{
			if (i == j)
				derivs[i][j] = 1;
		}
	}
	int g1, g2, g3, g4;
	g1 = 0; g2 = 0; g3 = 0; g4 = 0;
	for (int i = 0; i < size(derivs); i++)
	{
		if (double(i + 1) / size(derivs) <= 0.25)
		{
			SAD[g1].deriv = derivs[i]; g1++; 
		}
		else if (double(i + 1) / size(derivs) > 0.25 && double(i + 1) / size(derivs) <= 0.5)
		{
			sigmReverseAD[g2].deriv = derivs[i]; g2++; 
	    }
		else if (double(i + 1) / size(derivs) > 0.5 && double(i + 1) / size(derivs) <= 0.75)
		{
			rAD[g3].deriv = derivs[i]; g3++; 
		}
		else
		{
			tAD[g4].deriv = derivs[i]; g4++; 
		}
	}

	// compose 2 portfolios: one for double arguments and one for ForwardAD arguments
	vector<contract<double, double, double, double, double>> contractsDBL;
	vector<contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD>> contractsAD;
	for (int i = 0; i < a; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &sigma[i], &t[i], &r[i], T, K);
		contractsDBL.push_back(cDBL);
		contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD> cAD(&SAD[i], &sigmReverseAD[i], &tAD[i], &rAD[i], T, K);
		//contract<ForwardAD, double, double, double, ForwardAD> cAD(&SAD[i], &sigma[i], &t[i], &r[i], T, K);
		contractsAD.push_back(cAD);
	}

	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	portfolio<contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD>, ForwardAD> P_AD(contractsAD, "0");

	P_DBL.Price();
	P_AD.Price();

	cout << "Double Price: " << P_DBL.GetPrice() << endl;
	cout << "AD Price: " << P_AD.GetPrice() << endl;
	
	vector<double> aDeltas = P_DBL.Delta();
	vector<double> deltas = P_AD.GetDeltas();

	// delta, vega, rho, theta
	cout << "Delta: Analytical vs. AD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << aDeltas[i] << " - " << deltas[i] << " = " << aDeltas[i] - deltas[i] << endl;
	}
	cout << endl;

	vector<double> aVegas = P_DBL.Vega();
	vector<double> vegas = P_AD.GetVegas();

	cout << "Vega: Analytical vs. AD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << aVegas[i] << " - " << vegas[i] << " = " << aVegas[i] - vegas[i] << endl;
	}
	cout << endl;

	vector<double> aRhos = P_DBL.Rho();
	vector<double> rhos = P_AD.GetRhos();

	cout << "Rho: Analytical vs. AD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << aRhos[i] << " - " << rhos[i] << " = " << aRhos[i] - rhos[i] << endl;
	}
	cout << endl;

	vector<double> aThetas = P_DBL.Theta();
	vector<double> thetas = P_AD.GetThetas();

	cout << "Theta: Analytical vs. AD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << aThetas[i] << " - " << thetas[i] << " = " << aThetas[i] - thetas[i] << endl;
	}
	cout << endl;

	// delta, vega, rho, theta
	vector<double> fDeltas = P_DBL.FiniteDiff(S, 1);
	cout << "Delta: Analytical vs. FD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << aDeltas[i] << " - " << fDeltas[i] << " = " << aDeltas[i] - fDeltas[i] << endl;
	}
	cout << endl;

	vector<double> fVegas = P_DBL.FiniteDiff(sigma, 1);
	cout << "Vega: Analytical vs. FD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << aVegas[i] << " - " << fVegas[i] << " = " << aVegas[i] - fVegas[i] << endl;
	}
	cout << endl;

	vector<double> fRhos = P_DBL.FiniteDiff(r, 1);
	cout << "Rho: Analytical vs. FD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << aRhos[i] << " - " << fRhos[i] << " = " << aRhos[i] - fRhos[i] << endl;
	}
	cout << endl;

	vector<double> fThetas = P_DBL.FiniteDiff(t, 1);
	cout << "Theta: Analytical vs. FD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << aThetas[i] << " - " << fThetas[i] << " = " << aThetas[i] - fThetas[i] << endl;
	}

	clock_t time = clock();
	for (int i = 0; i < it; i++)
	{
		P_AD.Price();
	}
	cout << "AD calculation time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;

	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBL.Price();
	}
	cout << "DBL Price calculation time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;

	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBL.FiniteDiff(S, 1);
		P_DBL.FiniteDiff(sigma, 1);
		P_DBL.FiniteDiff(r, 1);
		P_DBL.FiniteDiff(t, 1);
	}
	cout << "FD calculation time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;

}

// contract <S, sigma, t, r, T, K>
void MonteCarloTest(const int it, const int a, const int n, const int m)
{
	// a - contracts
	// n - days
	// m - trajectories
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, a);
	vector<ForwardAD> SAD = VarVectorAD(aS, bS, a); // will not work if we have other AD variables except S
	double K = 300;
	double Sigma = 0.25;
	double r = 0.03;
	double T = 1;
	double t = 1. / 365;

	// compose 2 portfolios: one for double arguments and one for ForwardAD arguments
	vector<contract<double, double, double, double, double> > contractsDBL;
	vector<contract<ForwardAD, double, double, double, ForwardAD> > contractsAD;
	for (int i = 0; i < a; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<ForwardAD, double, double, double, ForwardAD> cAD(&SAD[i], &Sigma, &t, &r, T, K);
		contractsAD.push_back(cAD);
	}

	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	clock_t time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBL.Price();
	}
	cout << "DoubleAnalyticalPrice calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "DoubleAnalyticalPrice: " << P_DBL.GetPrice() << endl;
	cout << "AnalyticalDelta: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Delta()[i] << endl;
	}
	cout << endl;

	portfolio<contract<ForwardAD, double, double, double, ForwardAD>, ForwardAD> P_ADmc1(contractsAD, "MC1", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_ADmc1.Price();
	}
	cout << "AD MC1 calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "AD MC1 Price: " << P_ADmc1.GetPrice() << endl;
	cout << "AD MC1 Delta: " << endl;
	vector<double> deltas = P_ADmc1.GetDeltas();
	for (int i = 0; i < a; i++)
	{
		cout << deltas[i] << endl;
	}
	cout << endl;

	portfolio<contract<ForwardAD, double, double, double, ForwardAD>, ForwardAD> P_ADmc2(contractsAD, "MC2", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_ADmc2.Price();
	}
	cout << "AD MC2 calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "AD MC2 Price: " << P_ADmc2.GetPrice() << endl;
	cout << "AD MC2 Delta: " << endl;
	deltas = P_ADmc2.GetDeltas();
	for (int i = 0; i < a; i++)
	{
		cout << deltas[i] << endl;
	}
	cout << endl;

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc1(contractsDBL, "MC1", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBLmc1.Price();
	}
	cout << "Double MC1 Price calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double MC1 Price: " << P_DBLmc1.GetPrice() << endl;
	time = clock();
	deltas = P_DBLmc1.FiniteDiff(S, 1);
	cout << "Double MC1 FiniteDiff calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl << endl;
	cout << "Double MC1 FiniteDiff: " << endl;
	for (int i = 0; i < size(contractsDBL); i++)
	{
		cout << deltas[i] << endl;
	}

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc2(contractsDBL, "MC2", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBLmc2.Price();
	}
	cout << "Double MC2 Price calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double MC2 Price: " << P_DBLmc2.GetPrice() << endl;
	time = clock();
	deltas = P_DBLmc2.FiniteDiff(S, 1);
	cout << "Double MC2 FiniteDiff calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double MC2 FiniteDiff: " << endl;
	for (int i = 0; i < size(contractsDBL); i++)
	{
		cout << deltas[i] << endl;
	}
	
}

void ReverseADSpotTest(const int it, const int n)
{
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, n);
	vector<ReverseAD> SReverseAD = VarVectorReverseAD(aS, bS, n);
	double K = 300;
	double Sigma = 0.25;
	double r = 0.03;
	double T = 1;
	double t = 1. / 365;

	clock_t time = clock();

	// compose 2 portfolios: one for double arguments and one for ForwardAD arguments
	vector<contract<double, double, double, double, double> > contractsDBL;
	vector<contract<ReverseAD, double, double, double, ReverseAD> > contractsReverseAD;
	for (int i = 0; i < n; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<ReverseAD, double, double, double, ReverseAD> cReverseAD(&SReverseAD[i], &Sigma, &t, &r, T, K);
		contractsReverseAD.push_back(cReverseAD);
	}
	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	portfolio<contract<ReverseAD, double, double, double, ReverseAD>, ReverseAD> P_ReverseAD(contractsReverseAD, "0");

	cout << "DBL function calculation time = ";
	time = clock();
	for (int i = 0; i < it; ++i)
		P_DBL.Price();
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

	cout << "ReverseAD function calculation time = ";
	time = clock();
	for (int i = 0; i < it; ++i)
	{
		//for (int j = 0; j < n; ++j)
		//	SReverseAD[j].Reset();
		//priceAD.Reset();
		P_ReverseAD.Price();
	}
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

	// finite-difference derivatives
	vector<double> derivsFD(n);
	cout << "FD derivatives calculation time = ";
	time = clock();
	for (int i = 0; i < it; ++i)
		derivsFD = P_DBL.FiniteDiff(S, 1);
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl << endl;

	cout << "DBL function value = " << P_DBL.GetPrice() << endl;
	cout << "ReverseAD function value = " << P_ReverseAD.GetPrice() << endl;


	// printing out derivative values: AD versus analytical vs FD
	double totalDelta = 0;
	//cout << "Analytical deltas:" << endl;
	for (int i = 0; i < size(contractsDBL); ++i)
	{
		totalDelta += contractsDBL[i].Delta();
		//cout << contractsDBL[i].Delta() << "\t";
	}
	cout << endl << "Portfolio analytical total delta = " << totalDelta << endl;

	totalDelta = 0;
	//priceAD.SetGradient(1.0);
	//cout << "ReverseAD deltas:" << endl;
	vector<double> deltas;
	deltas = P_ReverseAD.GetDeltas();
	time = clock();
	for (int i = 0; i < size(contractsReverseAD); ++i) // in the case of this test we have n independent variables
	{
		//double tmp = SReverseAD[i].GetGradient();
		//totalDelta += tmp;
		totalDelta += deltas[i];
		//cout << tmp << "\t";
	}
	cout << endl << "Portfolio ReverseAD total delta = " << totalDelta << endl;
	//cout << "ReverseAD derivatives calculation time = ";
	//cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

	totalDelta = 0;
	//cout << "FD deltas:" << endl;
	for (int i = 0; i < n; ++i) // in the case of this test we have n independent variables
	{
		totalDelta += derivsFD[i];
		//cout << derivsFD[i] << "\t";
	}
	cout << endl << "Portfolio FD total delta = " << totalDelta << endl;


}

void ReverseADMonteCarloTest(const int it, const int a, const int n, const int m)
{
	// a - contracts
	// n - days
	// m - trajectories
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, a);
	vector<ReverseAD> SAD1 = VarVectorReverseAD(aS, bS, a); // will not work if we have other AD variables except S
	vector<ReverseAD> SAD2 = VarVectorReverseAD(aS, bS, a); // will not work if we have other AD variables except S
	double K = 300;
	double Sigma = 0.25;
	double r = 0.03;
	double T = 1;
	double t = 1. / 365;

	// compose 2 portfolios: one for double arguments and one for ForwardAD arguments
	vector<contract<double, double, double, double, double> > contractsDBL1, contractsDBL2;
	vector<contract<ReverseAD, double, double, double, ReverseAD> > contractsReverseAD1, contractsReverseAD2;
	for (int i = 0; i < a; i++)
	{
		contract<double, double, double, double, double> cDBL1(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL1.push_back(cDBL1);
		contract<double, double, double, double, double> cDBL2(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL2.push_back(cDBL2);
		contract<ReverseAD, double, double, double, ReverseAD> cAD1(&SAD1[i], &Sigma, &t, &r, T, K);
		contractsReverseAD1.push_back(cAD1);
		contract<ReverseAD, double, double, double, ReverseAD> cAD2(&SAD2[i], &Sigma, &t, &r, T, K);
		contractsReverseAD2.push_back(cAD2);
	}

	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL1, "0");
	clock_t time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBL.Price();
	}
	cout << "Double Analytical Price calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double Analytical Price: " << P_DBL.GetPrice() << endl;
	//cout << "Analytical Deltas: " << endl;
	double totalDelta = 0.0;
	for (int i = 0; i < a; i++)
	{
		double tmp = P_DBL.Delta()[i];
		totalDelta += tmp;
		//cout << tmp << endl;
	}
	//cout << endl;
	cout << "Cumulative Analytical Delta: " << totalDelta << endl;

	portfolio<contract<ReverseAD, double, double, double, ReverseAD>, ReverseAD> P_ADmc1(contractsReverseAD1, "MC1", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_ADmc1.Price(); // it is required to reset variables!
	}
	cout << "ReverseAD MC1 calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "ReverseAD MC1 Price: " << P_ADmc1.GetPrice() << endl;
	//cout << "ReverseAD MC1 Deltas: " << endl;
	//price1.SetGradient(1.0);
	vector<double> deltas = P_ADmc1.GetDeltas();
	totalDelta = 0.0;
	for (int i = 0; i < size(deltas); i++)
	{
		totalDelta += deltas[i];
		//cout << tmp << "\t";
	}
	//cout << endl;
	cout << "Cumulative ReverseAD MC1 Delta: " << totalDelta << endl;

	portfolio<contract<ReverseAD, double, double, double, ReverseAD>, ReverseAD> P_ADmc2(contractsReverseAD2, "MC2", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_ADmc2.Price();
	}
	cout << "ReverseAD MC2 calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "ReverseAD MC2 Price: " << P_ADmc2.GetPrice() << endl;
	//cout << "ReverseAD MC2 Deltas: " << endl;
	//price2.SetGradient(1.0);
	deltas = P_ADmc2.GetDeltas();
	totalDelta = 0.0;
	for (int i = 0; i < size(deltas); i++)
	{
		totalDelta += deltas[i];
		//cout << tmp << "\t";
	}
	//cout << endl;
	cout << "Cumulative ReverseAD MC2 Delta: " << totalDelta << endl;

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc1(contractsDBL1, "MC1", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBLmc1.Price();
	}
	cout << "Double MC1 Price calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double MC1 Price: " << P_DBLmc1.GetPrice() << endl;
	//cout << "Double MC1 FiniteDiff Deltas: " << endl;
	time = clock();
	deltas = P_DBLmc1.FiniteDiff(S, 1);
	cout << "Double MC1 FiniteDiff calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	totalDelta = 0.0;
	for (int i = 0; i < size(S); i++)
	{
		totalDelta += deltas[i];
		//cout << tmp << "\t";
	}
	cout << "Cumulative MC1 FiniteDiff Delta: " << totalDelta << endl;

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc2(contractsDBL2, "MC2", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBLmc2.Price();
	}
	cout << "Double MC2 Price calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double MC2 Price: " << P_DBLmc2.GetPrice() << endl;
	//cout << "Double MC2 FiniteDiff: " << endl;
	totalDelta = 0.0;
	time = clock();
	deltas = P_DBLmc2.FiniteDiff(S, 1);
	cout << "Double MC2 FiniteDiff calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	for (int i = 0; i < size(S); i++)
	{
		totalDelta += deltas[i];
		//cout << tmp << "\t";
	}
	cout << "Cumulative MC2 FiniteDiff Delta: " << totalDelta << endl;

}


