#pragma once
#include <RND.h>

void SpotTest(const int it, const int n)
{
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, n);
	vector<autodiff> SAD = VarVectorAD(aS, bS, n); // will not work if we have other AD variables except S
	double K = 300;
	double Sigma = 0.25;
	double r = 0.03;
	double T = 1;
	double t = 1. / 365;
	double priceDBL;
	autodiff priceAD;

	clock_t time = clock();

	// compose 2 portfolios: one for double arguments and one for autodiff arguments
	vector<contract<double, double, double, double, double> > contractsDBL;
	vector<contract<autodiff, double, double, double, autodiff> > contractsAD;
	for (int i = 0; i < n; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<autodiff, double, double, double, autodiff> cAD(&SAD[i], &Sigma, &t, &r, T, K);
		contractsAD.push_back(cAD);
	}
	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	portfolio<contract<autodiff, double, double, double, autodiff>, autodiff> P_AD(contractsAD, "0");
	cout << "AD function calculation time: ";
	time = clock();
	for (int i = 0; i < it; ++i)
		priceAD = P_AD.Price();
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "DBL function calculation time: ";
	time = clock();
	for (int i = 0; i < it; ++i)
		priceDBL = P_DBL.Price();
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

	// finite-difference derivatives
	vector<double> derivsFD(n);
	cout << "FD derivatives calculation time: ";
	time = clock();
	for (int i = 0; i < it; ++i)
		derivsFD = P_DBL.FiniteDiff(priceDBL, S, 1);
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

	cout << "DBL function value = " << priceDBL << endl;
	cout << "AD function value = " << priceAD.value << endl;


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
	for (int i = 0; i < n; ++i) // in the case of this test we have n independent variables
	{
		totalDelta += priceAD.deriv[i];
		cout << priceAD.deriv[i] << "\t";
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

}

void AllVarTest(const int it, const int a)
{
	double aS = 100; double bS = 500;
	double aSigma = 0.15; double bSigma = 0.3;
	double aR = 0.01; double bR = 0.04;
	double aT = 1. / 365;
	double K = 300;
	double T = 1;
	vector<double> S = VarVector(aS, bS, a); vector<autodiff>  SAD = VarVectorAD(aS, bS, a);
	vector<double> sigma = VarVector(aSigma, bSigma, a); vector<autodiff> sigmaAD = VarVectorAD(aSigma, bSigma, a);
	vector<double> 	r = VarVector(aR, bR, a); vector<autodiff>  rAD = VarVectorAD(aR, bR, a);
	vector<double> t = VarVector(aT, aT, a); vector<autodiff> tAD = VarVectorAD(aT, aT, a);
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
			sigmaAD[g2].deriv = derivs[i]; g2++; 
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

	// compose 2 portfolios: one for double arguments and one for autodiff arguments
	vector<contract<double, double, double, double, double>> contractsDBL;
	vector<contract<autodiff, autodiff, autodiff, autodiff, autodiff>> contractsAD;
	for (int i = 0; i < a; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &sigma[i], &t[i], &r[i], T, K);
		contractsDBL.push_back(cDBL);
		contract<autodiff, autodiff, autodiff, autodiff, autodiff> cAD(&SAD[i], &sigmaAD[i], &tAD[i], &rAD[i], T, K);
		//contract<autodiff, double, double, double, autodiff> cAD(&SAD[i], &sigma[i], &t[i], &r[i], T, K);
		contractsAD.push_back(cAD);
	}

	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	portfolio<contract<autodiff, autodiff, autodiff, autodiff, autodiff>, autodiff> P_AD(contractsAD, "0");

	cout << "Double Price: " << P_DBL.Price() << endl;
	cout << "AD Price: " << P_AD.Price().value << endl;
	

	// delta, vega, rho, theta
	cout << "Delta: Analytical vs. AD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Delta()[i] << " - " << P_AD.Price().deriv[i] << " = " << P_DBL.Delta()[i] - P_AD.Price().deriv[i] << endl;
	}
	cout << endl;
	cout << "Vega: Analytical vs. AD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Vega()[i] << " - " << P_AD.Price().deriv[i + a] << " = " << P_DBL.Vega()[i] - P_AD.Price().deriv[i + a] << endl;
	}
	cout << endl;
	cout << "Rho: Analytical vs. AD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Rho()[i] << " - " << P_AD.Price().deriv[i + 2*a] << " = " << P_DBL.Rho()[i] - P_AD.Price().deriv[i + 2*a] << endl;
	}
	cout << endl;
	cout << "Theta: Analytical vs. AD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Theta()[i] << " - " << P_AD.Price().deriv[i + 3*a] << " = " << P_DBL.Theta()[i] - P_AD.Price().deriv[i + 3*a] << endl;
	}
	cout << endl;
	// delta, vega, rho, theta
	cout << "Delta: Analytical vs. FD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Delta()[i] << " - " << P_DBL.FiniteDiff(P_DBL.Price(), S, 1)[i] << " = " << P_DBL.Delta()[i] - P_DBL.FiniteDiff(P_DBL.Price(), S, 1)[i] << endl;
	}
	cout << endl;
	cout << "Vega: Analytical vs. FD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Vega()[i] << " - " << P_DBL.FiniteDiff(P_DBL.Price(), sigma, 1)[i] << " = " << P_DBL.Vega()[i] - P_DBL.FiniteDiff(P_DBL.Price(), sigma, 1)[i] << endl;
	}
	cout << endl;
	cout << "Rho: Analytical vs. FD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Rho()[i] << " - " << P_DBL.FiniteDiff(P_DBL.Price(), r, 1)[i] << " = " << P_DBL.Rho()[i] - P_DBL.FiniteDiff(P_DBL.Price(), r, 1)[i] << endl;
	}
	cout << endl;
	cout << "Theta: Analytical vs. FD: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Theta()[i] << " - " << P_DBL.FiniteDiff(P_DBL.Price(), t, 1)[i] << " = " << P_DBL.Theta()[i] - P_DBL.FiniteDiff(P_DBL.Price(), t, 1)[i] << endl;
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
		P_DBL.FiniteDiff(P_DBL.Price(), S, 1);
		P_DBL.FiniteDiff(P_DBL.Price(), sigma, 1);
		P_DBL.FiniteDiff(P_DBL.Price(), r, 1);
		P_DBL.FiniteDiff(P_DBL.Price(), t, 1);
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
	vector<autodiff> SAD = VarVectorAD(aS, bS, a); // will not work if we have other AD variables except S
	double K = 300;
	double Sigma = 0.25;
	double r = 0.03;
	double T = 1;
	double t = 1. / 365;

	// compose 2 portfolios: one for double arguments and one for autodiff arguments
	vector<contract<double, double, double, double, double> > contractsDBL;
	vector<contract<autodiff, double, double, double, autodiff> > contractsAD;
	for (int i = 0; i < a; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<autodiff, double, double, double, autodiff> cAD(&SAD[i], &Sigma, &t, &r, T, K);
		contractsAD.push_back(cAD);
	}

	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	clock_t time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBL.Price();
	}
	cout << "DoubleAnalyticalPrice calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "DoubleAnalyticalPrice: " << P_DBL.Price() << endl;
	cout << "AnalyticalDelta: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Delta()[i] << endl;
	}
	cout << endl;

	portfolio<contract<autodiff, double, double, double, autodiff>, autodiff> P_ADmc1(contractsAD, "MC1", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_ADmc1.Price();
	}
	cout << "AD MC1 calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "AD MC1 Price: " << P_ADmc1.Price().value << endl;
	cout << "AD MC1 Delta: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_ADmc1.Price().deriv[i] << endl;
	}
	cout << endl;

	portfolio<contract<autodiff, double, double, double, autodiff>, autodiff> P_ADmc2(contractsAD, "MC2", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_ADmc1.Price();
	}
	cout << "AD MC2 calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "AD MC2 Price: " << P_ADmc2.Price().value << endl;
	cout << "AD MC2 Delta: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_ADmc2.Price().deriv[i] << endl;
	}
	cout << endl;

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc1(contractsDBL, "MC1", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBLmc1.Price();
	}
	cout << "Double MC1 Price calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double MC1 Price: " << P_DBLmc1.Price() << endl;
	cout << "Double MC1 FiniteDiff: " << endl;
	time = clock();
	for (int i = 0; i < it; i++)
	{
		cout << P_DBLmc1.FiniteDiff(P_DBLmc1.Price(), S, 1)[i] << endl;
	}
	cout << "Double MC1 FiniteDiff calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl << endl;

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc2(contractsDBL, "MC2", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBLmc2.Price();
	}
	cout << "Double MC2 Price calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double MC2 Price: " << P_DBLmc2.Price() << endl;
	cout << "Double MC2 FiniteDiff: " << endl;
	time = clock();
	for (int i = 0; i < it; i++)
	{
		cout << P_DBLmc2.FiniteDiff(P_DBLmc2.Price(), S, 1)[i] << endl;
	}
	cout << "Double MC2 FiniteDiff calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

}

void AADSpotTest(const int it, const int n)
{
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, n);
	vector<aad> SAAD = VarVectorAAD(aS, bS, n);
	double K = 300;
	double Sigma = 0.25;
	double r = 0.03;
	double T = 1;
	double t = 1. / 365;
	double priceDBL;
	aad priceAD;

	clock_t time = clock();

	// compose 2 portfolios: one for double arguments and one for autodiff arguments
	vector<contract<double, double, double, double, double> > contractsDBL;
	vector<contract<aad, double, double, double, aad> > contractsAD;
	for (int i = 0; i < n; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<aad, double, double, double, aad> cAD(&SAAD[i], &Sigma, &t, &r, T, K);
		contractsAD.push_back(cAD);
	}
	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	portfolio<contract<aad, double, double, double, aad>, aad> P_AD(contractsAD, "0");

	cout << "DBL function calculation time" << endl;
	time = clock();
	for (int i = 0; i < it; ++i)
		priceDBL = P_DBL.Price();
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl << endl;

	cout << "AAD function calculation time" << endl;
	time = clock();
	for (int i = 0; i < it; ++i)
	{
		for (int j = 0; j < n; ++j)
			SAAD[j].Reset();
		priceAD.Reset();
		priceAD = P_AD.Price();
	}
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl << endl;

	// finite-difference derivatives
	vector<double> derivsFD(n);
	cout << "FD derivatives calculation time" << endl;
	time = clock();
	for (int i = 0; i < it; ++i)
		derivsFD = P_DBL.FiniteDiff(priceDBL, S, 1);
	cout << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

	cout << "DBL function value = " << priceDBL << endl;
	cout << "AAD function value = " << priceAD.GetVar()->GetValue() << endl;


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
	//priceAD.SetGradient(1.0);
	cout << "AAD deltas:" << endl;
	for (int i = 0; i < n; ++i) // in the case of this test we have n independent variables
	{
		double tmp = SAAD[i].GetGradient();
		totalDelta += tmp;
		cout << tmp << "\t";
	}
	cout << endl << "Portfolio AAD total delta = " << totalDelta << endl;

	totalDelta = 0;
	cout << "FD deltas:" << endl;
	for (int i = 0; i < n; ++i) // in the case of this test we have n independent variables
	{
		totalDelta += derivsFD[i];
		cout << derivsFD[i] << "\t";
	}
	cout << endl << "Portfolio FD total delta = " << totalDelta << endl;

	system("pause");
}

void AADMonteCarloTest(const int it, const int a, const int n, const int m)
{
	// a - contracts
	// n - days
	// m - trajectories
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, a);
	vector<aad> SAD = VarVectorAAD(aS, bS, a); // will not work if we have other AD variables except S
	double K = 300;
	double Sigma = 0.25;
	double r = 0.03;
	double T = 1;
	double t = 1. / 365;
	aad price1, price2;

	// compose 2 portfolios: one for double arguments and one for autodiff arguments
	vector<contract<double, double, double, double, double> > contractsDBL;
	vector<contract<aad, double, double, double, aad> > contractsAD;
	for (int i = 0; i < a; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<aad, double, double, double, aad> cAD(&SAD[i], &Sigma, &t, &r, T, K);
		contractsAD.push_back(cAD);
	}

	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	clock_t time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBL.Price();
	}
	cout << "DoubleAnalyticalPrice calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "DoubleAnalyticalPrice: " << P_DBL.Price() << endl;
	cout << "AnalyticalDelta: " << endl;
	for (int i = 0; i < a; i++)
	{
		cout << P_DBL.Delta()[i] << endl;
	}
	cout << endl;

	portfolio<contract<aad, double, double, double, aad>, aad> P_ADmc1(contractsAD, "MC1", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		price1 = P_ADmc1.Price();
	}
	cout << "AD MC1 calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "AD MC1 Price: " << price1.GetVar()->GetValue() << endl;
	cout << "AD MC1 Delta: " << endl;
	//price1.SetGradient(1.0);
	double totalDelta = 0.0;
	for (int i = 0; i < a; i++)
	{
		double tmp = SAD[i].GetGradient();
		totalDelta += tmp;
		cout << tmp << "\t";
	}
	cout << endl;

	portfolio<contract<aad, double, double, double, aad>, aad> P_ADmc2(contractsAD, "MC2", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		price2 = P_ADmc1.Price();
	}
	cout << "AD MC2 calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "AD MC2 Price: " << price2.GetVar()->GetValue() << endl;
	cout << "AD MC2 Delta: " << endl;
	//price2.SetGradient(1.0);
	totalDelta = 0.0;
	for (int i = 0; i < a; i++)
	{
		double tmp = SAD[i].GetGradient();
		totalDelta += tmp;
		cout << tmp << "\t";
	}
	cout << endl;

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc1(contractsDBL, "MC1", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBLmc1.Price();
	}
	cout << "Double MC1 Price calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double MC1 Price: " << P_DBLmc1.Price() << endl;
	cout << "Double MC1 FiniteDiff: " << endl;
	time = clock();
	for (int i = 0; i < it; i++)
	{
		cout << P_DBLmc1.FiniteDiff(P_DBLmc1.Price(), S, 1)[i] << endl;
	}
	cout << "Double MC1 FiniteDiff calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl << endl;

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc2(contractsDBL, "MC2", n, m);
	time = clock();
	for (int i = 0; i < it; i++)
	{
		P_DBLmc2.Price();
	}
	cout << "Double MC2 Price calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;
	cout << "Double MC2 Price: " << P_DBLmc2.Price() << endl;
	cout << "Double MC2 FiniteDiff: " << endl;
	time = clock();
	for (int i = 0; i < it; i++)
	{
		cout << P_DBLmc2.FiniteDiff(P_DBLmc2.Price(), S, 1)[i] << endl;
	}
	cout << "Double MC2 FiniteDiff calculation time: " << (clock() - time) / (double)CLOCKS_PER_SEC << endl;

}


