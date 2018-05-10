#pragma once
#include <RND.h>


void ReverseTest(const int it, const int a, const int n, const int m)
{
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, a);
	vector<ReverseAD> SRAD = VarVectorReverseAD(aS, bS, a);
	double K = 300;
	double Sigma = 0.25;
	vector<double> derivsig(a);
	derivsig.push_back(1);
	ReverseAD SigmaRAD(Sigma);
	ReverseAD SigmaRAD1(Sigma);
	ReverseAD SigmaRAD2(Sigma);
	double r = 0.03;
	vector<double> derivr(a + 2);
	derivr.push_back(1);

	ReverseAD rRAD(0.03);
	ReverseAD rRAD1(0.03);
	ReverseAD rRAD2(0.03);
	double T = 1;
	double t = 1. / 365;
	vector<double> derivt(a + 1);
	derivt.push_back(1);

	ReverseAD tRAD(1. / 365);
	ReverseAD tRAD1(1. / 365);
	ReverseAD tRAD2(1. / 365);


	vector<contract<double, double, double, double, double>> contractsDBL, contractsDBL1, contractsDBL2;
	vector<contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD>> contractsReverseAD, contractsReverseAD1, contractsReverseAD2;

	for (int i = 0; i < a; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<double, double, double, double, double> cDBL1(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL1.push_back(cDBL1);
		contract<double, double, double, double, double> cDBL2(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL2.push_back(cDBL2);

		contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD> cRAD(&SRAD[i], &SigmaRAD, &tRAD, &rRAD, T, K);
		contractsReverseAD.push_back(cRAD);
		contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD> cRAD1(&SRAD[i], &SigmaRAD1, &tRAD1, &rRAD1, T, K);
		contractsReverseAD1.push_back(cRAD1);
		contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD> cRAD2(&SRAD[i], &SigmaRAD2, &tRAD2, &rRAD2, T, K);
		contractsReverseAD2.push_back(cRAD2);
	}
	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	portfolio<contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD>, ReverseAD> P_RAD(contractsReverseAD, "0");

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc1(contractsDBL1, "MC1", n, m);
	portfolio<contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD>, ReverseAD> P_RADmc1(contractsReverseAD1, "MC1", n, n);

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc2(contractsDBL2, "MC2", n, m);
	portfolio<contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD>, ReverseAD> P_RADmc2(contractsReverseAD2, "MC2", n, n);

	clock_t time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_RAD.Price();
	}
	cout << "ReverseAD 0 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cout << "ReverseAD 0 Price: " << P_RAD.GetPrice() << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_RADmc1.Price();
	}
	cout << "ReverseAD MC1 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cout << "ReverseAD MC1 Price: " << P_RADmc1.GetPrice() << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_RADmc2.Price();
	}
	cout << "ReverseAD MC2 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cout << "ReverseAD MC2 Price: " << P_RADmc2.GetPrice() << endl;

	P_RAD.Price();
	P_RADmc1.Price();

	P_RADmc2.Price();


	time = clock();
	for (int i = 0; i < 1000 * it; ++i)
	{
		P_DBL.Price();
	}
	cout << "DBL Price 0 time: " << double(clock() - time) / CLOCKS_PER_SEC / 1000 << endl;
	cout << "Theory Price: " << P_DBL.GetPrice() << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_DBLmc1.Price();
	}
	cout << "DBL Price MC1 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_DBLmc2.Price();
	}
	cout << "DBL Price MC2 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;



	double totalDelta = 0;	double totalVega = 0;	double totalTheta = 0;	double totalRho = 0;
	vector<double> deltas = P_DBL.Delta();
	vector<double> vegas = P_DBL.Vega();
	vector<double> thetas = P_DBL.Theta();
	vector<double> rhos = P_DBL.Rho();
	for (int i = 0; i < a; ++i)
	{
		totalDelta += deltas[i];
		totalVega += vegas[i];
		totalTheta += thetas[i];
		totalRho += rhos[i];
	}
	cout << "Theory total delta = " << totalDelta << endl;
	cout << "Theory total vega = " << totalVega << endl;
	cout << "Theory total theta = " << totalTheta << endl;
	cout << "Theory total rho = " << totalRho << endl;


	
	double totalDeltaR = 0;	double totalVegaR = 0;	double totalThetaR = 0;	double totalRhoR = 0;
	vector<double> deltasR = P_RAD.GetDeltas();
	vector<double> vegasR = P_RAD.GetVegas();
	vector<double> thetasR = P_RAD.GetThetas();
	vector<double> rhosR = P_RAD.GetRhos();
	for (int i = 0; i < a; ++i)
	{
		totalDeltaR += deltasR[i];
		totalVegaR += vegasR[i];
		totalThetaR += thetasR[i];
		totalRhoR += rhosR[i];
	}
	cout << "Error ReverseAD 0 total delta = " << abs(totalDelta - totalDeltaR) << endl;
	cout << "ReverseAD 0 total vega = " << totalVegaR << endl;
	cout << "ReverseAD 0 total theta = " << totalThetaR << endl;
	cout << "ReverseAD 0 total rho = " << totalRhoR << endl;

	double totalDeltaR1 = 0;	double totalVegaR1 = 0;	double totalThetaR1 = 0;	double totalRhoR1 = 0;
	vector<double> deltasR1 = P_RADmc1.GetDeltas();
	vector<double> vegasR1 = P_RADmc1.GetVegas();
	vector<double> thetasR1 = P_RADmc1.GetThetas();
	vector<double> rhosR1 = P_RADmc1.GetRhos();
	for (int i = 0; i < a; ++i)
	{
		totalDeltaR1 += deltasR1[i];
		totalVegaR1 += vegasR1[i];
		totalThetaR1 += thetasR1[i];
		totalRhoR1 += rhosR1[i];
	}
	cout << "Error ReverseAD MC1 total delta = " << abs(totalDelta - totalDeltaR1) << endl;
	cout << "ReverseAD MC1 total vega = " << totalVegaR1 << endl;
	cout << "ReverseAD MC1 total theta = " << totalThetaR1 << endl;
	cout << "ReverseAD MC1 total rho = " << totalRhoR1 << endl;

	double totalDeltaR2 = 0;	double totalVegaR2 = 0;	double totalThetaR2 = 0;	double totalRhoR2 = 0;
	vector<double> deltasR2 = P_RADmc2.GetDeltas();
	vector<double> vegasR2 = P_RADmc2.GetVegas();
	vector<double> thetasR2 = P_RADmc2.GetThetas();
	vector<double> rhosR2 = P_RADmc2.GetRhos();
	for (int i = 0; i < a; ++i)
	{
		totalDeltaR2 += deltasR2[i];
		totalVegaR2 += vegasR2[i];
		totalThetaR2 += thetasR2[i];
		totalRhoR2 += rhosR2[i];
	}
	cout << "Error ReverseAD MC2 total delta = " << abs(totalDelta - totalDeltaR2) << endl;
	cout << "ReverseAD MC2 total vega = " << totalVegaR2 << endl;
	cout << "ReverseAD MC2 total theta = " << totalThetaR2 << endl;
	cout << "ReverseAD MC2 total rho = " << totalRhoR2 << endl;


	int bump = 1;
	vector<double>sgm(1, Sigma);
	vector<double>tt(1, t);
	vector<double>rr(1, r);
	double totalDeltafin = 0;
	vector<double> deltasfin = P_DBL.FiniteDiff(S, bump);
	for (int i = 0; i < a; ++i)
	{
		totalDeltafin += deltasfin[i];
	}
	cout << "Theory - FiniteDiff 0 total delta = " << abs(totalDelta - totalDeltafin) << endl;
	cout << "Theory - FiniteDiff 0 total vega = " << abs(totalVega - P_DBL.FiniteDiff(sgm, bump)[0]) << endl;
	cout << "Theory - FiniteDiff 0 total theta = " << abs(totalTheta - P_DBL.FiniteDiff(tt, bump)[0]) << endl;
	cout << "Theory - FiniteDiff 0 total rho = " << abs(totalRho - P_DBL.FiniteDiff(tt, bump)[0]) << endl;

	double totalDeltafin1 = 0;
	vector<double> deltasfin1 = P_DBLmc1.FiniteDiff(S, bump);
	for (int i = 0; i < a; ++i)
	{
		totalDeltafin1 += deltasfin1[i];
	}
	cout << "Theory - FiniteDiff MC1 total delta = " << abs(totalDelta - totalDeltafin1) << endl;
	cout << "Theory - FiniteDiff MC1 total vega = " << abs(totalVega - P_DBLmc1.FiniteDiff(sgm, bump)[0]) << endl;
	cout << "Theory - FiniteDiff MC1 total theta = " << abs(totalTheta - P_DBLmc1.FiniteDiff(tt, bump)[0]) << endl;
	cout << "Theory - FiniteDiff MC1 total rho = " << abs(totalRho - P_DBLmc1.FiniteDiff(tt, bump)[0]) << endl;

	double totalDeltafin2 = 0;
	vector<double> deltasfin2 = P_DBLmc1.FiniteDiff(S, bump);
	for (int i = 0; i < a; ++i)
	{
		totalDeltafin2 += deltasfin2[i];
	}
	cout << "Theory - FiniteDiff MC2 total delta = " << abs(totalDelta - totalDeltafin2) << endl;
	cout << "Theory - FiniteDiff MC2 total vega = " << abs(totalVega - P_DBLmc2.FiniteDiff(sgm, bump)[0]) << endl;
	cout << "Theory - FiniteDiff MC2 total theta = " << abs(totalTheta - P_DBLmc2.FiniteDiff(tt, bump)[0]) << endl;
	cout << "Theory - FiniteDiff MC2 total rho = " << abs(totalRho - P_DBLmc2.FiniteDiff(tt, bump)[0]) << endl;

	time = clock();
	for (int i = 0; i < 1000 * it; ++i)
	{
		P_DBL.FiniteDiff(S, bump);
		P_DBL.FiniteDiff(sgm, bump);
		P_DBL.FiniteDiff(tt, bump);
		P_DBL.FiniteDiff(rr, bump);
	}
	cout << "DBL FiniteDiff 0 time: " << double(clock() - time) / CLOCKS_PER_SEC / 1000 << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_DBLmc1.FiniteDiff(S, bump);
		P_DBLmc1.FiniteDiff(sgm, bump);
		P_DBLmc1.FiniteDiff(tt, bump);
		P_DBLmc1.FiniteDiff(rr, bump);
	}
	cout << "DBL FiniteDiff MC1 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_DBLmc2.FiniteDiff(S, bump);
		P_DBLmc2.FiniteDiff(sgm, bump);
		P_DBLmc2.FiniteDiff(tt, bump);
		P_DBLmc2.FiniteDiff(rr, bump);
	}
	cout << "DBL FiniteDiff MC2 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;

	//cout << "------------0-------------" << endl;
	//cout << "DBL function value = " << P_DBL.GetPrice() << endl;
	//cout << "ReverseAD function value = " << P_RAD.GetPrice() << endl;
	//cout << "ForwardAD function value = " << P_FAD.GetPrice() << endl;
	//cout << "DBL function value - ReverseAD function value = " << abs(P_DBL.GetPrice()- P_RAD.GetPrice()) << endl;
	//cout << "DBL function value - ForwardAD function value = " << abs(P_DBL.GetPrice() - P_FAD.GetPrice()) << endl;
	//cout << "------------MC1-------------" << endl;
	//cout << "DBL function value = " << P_DBLmc1.GetPrice() << endl;
	//cout << "Theory - DBL function value = " << P_DBL.GetPrice() - P_DBLmc1.GetPrice() << endl;
	//cout << "ReverseAD function value = " << P_RADmc1.GetPrice() << endl;
	//cout << "ForwardAD function value = " << P_FADmc1.GetPrice() << endl;
	//cout << "Theory - ReverseAD function value = " << abs(P_DBL.GetPrice() - P_RADmc1.GetPrice()) << endl;
	//cout << "Theory - ForwardAD function value = " << abs(P_DBL.GetPrice() - P_FADmc1.GetPrice()) << endl;
	//cout << "------------MC2-------------" << endl;
	//cout << "DBL function value = " << P_DBLmc2.GetPrice() << endl;
	//cout << "Theory - DBL function value = " << P_DBL.GetPrice() - P_DBLmc2.GetPrice() << endl;
	//cout << "ReverseAD function value = " << P_RADmc2.GetPrice() << endl;
	//cout << "ForwardAD function value = " << P_FADmc2.GetPrice() << endl;
	//cout << "Theory - ReverseAD function value = " << abs(P_DBL.GetPrice() - P_RADmc2.GetPrice()) << endl;
	//cout << "Theory - ForwardAD function value = " << abs(P_DBL.GetPrice() - P_FADmc2.GetPrice()) << endl;
}
void dissertation(const int it, const int a, const int n, const int m)
{
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, a);
	vector<ForwardAD> SFAD = VarVectorForwardAD(aS, bS, a);
	vector<ReverseAD> SRAD = VarVectorReverseAD(aS, bS, a);
	double K = 300;
	double Sigma = 0.25; 
	vector<double> derivsig(a);
	derivsig.push_back(1);
	ForwardAD SigmaFAD(Sigma, derivsig);
	ReverseAD SigmaRAD(0.25);
	ReverseAD SigmaRAD1(0.25);
	ReverseAD SigmaRAD2(0.25);
	double r = 0.03;
	vector<double> derivr(a + 2);
	derivr.push_back(1);
	ForwardAD rFAD(r, derivr);
	ReverseAD rRAD(0.03);
	ReverseAD rRAD1(0.03);
	ReverseAD rRAD2(0.03);
	double T = 1;
	double t = 1. / 365;
	vector<double> derivt(a + 1);
	derivt.push_back(1);
	ForwardAD tFAD(t, derivt);
	ReverseAD tRAD(1. / 365);
	ReverseAD tRAD1(1. / 365);
	ReverseAD tRAD2(1. / 365);
	for (int i = 0; i < a; i++)
	{
		for (int j = 1; j <= 3; j++)
		{
			SFAD[i].deriv.push_back(0);
		}
	}
	for (int j = 1; j <= 2; j++)
	{
		SigmaFAD.deriv.push_back(0);
	}
	tFAD.deriv.push_back(0);

	vector<contract<double, double, double, double, double>> contractsDBL, contractsDBL1, contractsDBL2;
	vector<contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD>> contractsForwardAD, contractsForwardAD1, contractsForwardAD2;
	vector<contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD>> contractsReverseAD, contractsReverseAD1, contractsReverseAD2;

	for (int i = 0; i < a; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<double, double, double, double, double> cDBL1(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL1.push_back(cDBL1);
		contract<double, double, double, double, double> cDBL2(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL2.push_back(cDBL2);
		contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD> cFAD(&SFAD[i], &SigmaFAD, &tFAD, &rFAD, T, K);
		contractsForwardAD.push_back(cFAD);
		contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD> cFAD1(&SFAD[i], &SigmaFAD, &tFAD, &rFAD, T, K);
		contractsForwardAD1.push_back(cFAD1);
		contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD> cFAD2(&SFAD[i], &SigmaFAD, &tFAD, &rFAD, T, K);
		contractsForwardAD2.push_back(cFAD2);
		contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD> cRAD(&SRAD[i], &SigmaRAD, &tRAD, &rRAD, T, K);
		contractsReverseAD.push_back(cRAD);
		contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD> cRAD1(&SRAD[i], &SigmaRAD1, &tRAD1, &rRAD1, T, K);
		contractsReverseAD1.push_back(cRAD1);
		contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD> cRAD2(&SRAD[i], &SigmaRAD2, &tRAD2, &rRAD2, T, K);
		contractsReverseAD2.push_back(cRAD2);
	}
	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
	portfolio<contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD>, ReverseAD> P_RAD(contractsReverseAD, "0");
	portfolio<contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD>, ForwardAD> P_FAD(contractsForwardAD, "0");

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc1(contractsDBL1, "MC1", n, m);
	portfolio<contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD>, ReverseAD> P_RADmc1(contractsReverseAD1, "MC1", n, n);
	portfolio<contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD>, ForwardAD> P_FADmc1(contractsForwardAD1, "MC1", n, m);

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc2(contractsDBL2, "MC2", n, m);
	portfolio<contract<ReverseAD, ReverseAD, ReverseAD, ReverseAD, ReverseAD>, ReverseAD> P_RADmc2(contractsReverseAD2, "MC2", n, n);
	portfolio<contract<ForwardAD, ForwardAD, ForwardAD, ForwardAD, ForwardAD>, ForwardAD> P_FADmc2(contractsForwardAD2, "MC2", n, m);

	clock_t time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_RAD.Price();
	}
	cout << "ReverseAD 0 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cout << "ReverseAD 0 Price: " << P_RAD.GetPrice() << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_RADmc1.Price();
	}
	cout << "ReverseAD MC1 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cout << "ReverseAD MC1 Price: " << P_RADmc1.GetPrice() << endl;
	
	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_RADmc2.Price();
	}
	cout << "ReverseAD MC2 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cout << "ReverseAD MC2 Price: " << P_RADmc2.GetPrice() << endl;


	//time = clock();
	//for (int i = 0; i < it; ++i)
	//{
	//	P_FAD.Price();
	//}
	//cout << "ForwardAD 0 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;
	//cout << "ForwardAD 0 Price: " << P_FAD.GetPrice() << endl;


	//time = clock();
	//for (int i = 0; i < it; ++i)
	//{
	//	P_FADmc1.Price();
	//}
	//cout << "ForwardAD MC1 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;
	//cout << "ForwardAD MC1 Price: " << P_FADmc1.GetPrice() << endl;


	//time = clock();
	//for (int i = 0; i < it; ++i)
	//{
	//	P_FADmc2.Price();
	//}
	//cout << "ForwardAD MC2 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;
	//cout << "ForwardAD MC2 Price: " << P_FADmc2.GetPrice() << endl;


	time = clock();
	for (int i = 0; i < 1000*it; ++i)
	{
		P_DBL.Price();
	}
	cout << "DBL Price 0 time: " << double(clock() - time) / CLOCKS_PER_SEC / 1000 << endl;
	cout << "Theory Price: " << P_DBL.GetPrice() << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_DBLmc1.Price();
	}
	cout << "DBL Price MC1 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_DBLmc2.Price();
	}
	cout << "DBL Price MC2 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;



	double totalDelta = 0;	double totalVega = 0;	double totalTheta = 0;	double totalRho = 0;
	vector<double> deltas = P_DBL.Delta();
	vector<double> vegas = P_DBL.Vega();
	vector<double> thetas = P_DBL.Theta();
	vector<double> rhos = P_DBL.Rho();
	for (int i = 0; i < a; ++i)
	{
		totalDelta += deltas[i];
		totalVega += vegas[i];
		totalTheta += thetas[i];
		totalRho += rhos[i];
	}
	cout  << "Theory total delta = " << totalDelta << endl;
	cout  << "Theory total vega = " << totalVega << endl;
	cout  << "Theory total theta = " << totalTheta << endl;
	cout  << "Theory total rho = " << totalRho << endl;


	//double totalDeltaF = 0;	double totalVegaF = 0;	double totalThetaF = 0;	double totalRhoF = 0;
	//vector<double> deltasF = P_FAD.GetDeltas();
	//vector<double> vegasF = P_FAD.GetVegas();
	//vector<double> thetasF = P_FAD.GetThetas();
	//vector<double> rhosF = P_FAD.GetRhos();
	//for (int i = 0; i < a; ++i)
	//{
	//	totalDeltaF += deltasF[i];
	//	totalVegaF += vegasF[i];
	//	totalThetaF += thetasF[i];
	//	totalRhoF += rhosF[i];
	//}
	//cout << "Error ForwardAD 0 total delta = " << abs(totalDelta - totalDeltaF) << endl;
	//cout << "Error ForwardAD 0 total vega = " << abs(totalVega - totalVega) << endl;
	//cout << "Error ForwardAD 0 total theta = " << abs(totalTheta - totalThetaF) << endl;
	//cout << "Error ForwardAD 0 total rho = " << abs(totalRho - totalRhoF) << endl;

	//double totalDeltaF1 = 0;	double totalVegaF1 = 0;	double totalThetaF1 = 0;	double totalRhoF1 = 0;
	//vector<double> deltasF1 = P_FADmc1.GetDeltas();
	//vector<double> vegasF1 = P_FADmc1.GetVegas();
	//vector<double> thetasF1 = P_FADmc1.GetThetas();
	//vector<double> rhosF1 = P_FADmc1.GetRhos();
	//for (int i = 0; i < a; ++i)
	//{
	//	totalDeltaF1 += deltasF1[i];
	//	totalVegaF1 += vegasF1[i];
	//	totalThetaF1 += thetasF1[i];
	//	totalRhoF1 += rhosF1[i];
	//}
	//cout << "Error ForwardAD MC1 total delta = " << abs(totalDelta - totalDeltaF1) << endl;
	//cout << "Error ForwardAD MC1 total vega = " << abs(totalVega - totalVegaF1) << endl;
	//cout << "Error ForwardAD MC1 total theta = " << abs(totalTheta - totalThetaF1) << endl;
	//cout << "Error ForwardAD MC1 total rho = " << abs(totalRho - totalRhoF1) << endl;

	//double totalDeltaF2 = 0;	double totalVegaF2 = 0;	double totalThetaF2 = 0;	double totalRhoF2 = 0;
	//vector<double> deltasF2 = P_FADmc2.GetDeltas();
	//vector<double> vegasF2 = P_FADmc2.GetVegas();
	//vector<double> thetasF2 = P_FADmc2.GetThetas();
	//vector<double> rhosF2 = P_FADmc2.GetRhos();
	//for (int i = 0; i < a; ++i)
	//{
	//	totalDeltaF2 += deltasF2[i];
	//	totalVegaF2 += vegasF2[i];
	//	totalThetaF2 += thetasF2[i];
	//	totalRhoF2 += rhosF2[i];
	//}
	//cout  << "Error ForwardAD MC2 total delta = " << abs(totalDelta - totalDeltaF2) << endl;
	//cout  << "Error ForwardAD MC2 total vega = " << abs(totalVega - totalVegaF2) << endl;
	//cout  << "Error ForwardAD MC2 total theta = " << abs(totalTheta - totalThetaF2) << endl;
	//cout  << "Error ForwardAD MC2 total rho = " << abs(totalRho - totalRhoF2) << endl;

	double totalDeltaR = 0;	double totalVegaR = 0;	double totalThetaR = 0;	double totalRhoR = 0;
	vector<double> deltasR = P_RAD.GetDeltas();
	vector<double> vegasR = P_RAD.GetVegas();
	vector<double> thetasR = P_RAD.GetThetas();
	vector<double> rhosR = P_RAD.GetRhos();
	for (int i = 0; i < a; ++i)
	{
		totalDeltaR += deltasR[i];
		totalVegaR += vegasR[i];
		totalThetaR += thetasR[i];
		totalRhoR += rhosR[i];
	}
	cout  << "Error ReverseAD 0 total delta = " << abs(totalDelta - totalDeltaR) << endl;
	cout  << "ReverseAD 0 total vega = " << totalVegaR << endl;
	cout  << "ReverseAD 0 total theta = " << totalThetaR << endl;
	cout  << "ReverseAD 0 total rho = " << totalRhoR << endl;

	double totalDeltaR1 = 0;	double totalVegaR1 = 0;	double totalThetaR1 = 0;	double totalRhoR1 = 0;
	vector<double> deltasR1 = P_RADmc1.GetDeltas();
	vector<double> vegasR1 = P_RADmc1.GetVegas();
	vector<double> thetasR1 = P_RADmc1.GetThetas();
	vector<double> rhosR1 = P_RADmc1.GetRhos();
	for (int i = 0; i < a; ++i)
	{
		totalDeltaR1 += deltasR1[i];
		totalVegaR1 += vegasR1[i];
		totalThetaR1 += thetasR1[i];
		totalRhoR1 += rhosR1[i];
	}
	cout  << "Error ReverseAD MC1 total delta = " << abs(totalDelta - totalDeltaR1) << endl;
	cout  << "ReverseAD MC1 total vega = " << totalVegaR1 << endl;
	cout  << "ReverseAD MC1 total theta = " << totalThetaR1 << endl;
	cout  << "ReverseAD MC1 total rho = " << totalRhoR1 << endl;

	double totalDeltaR2 = 0;	double totalVegaR2 = 0;	double totalThetaR2 = 0;	double totalRhoR2 = 0;
	vector<double> deltasR2 = P_RADmc2.GetDeltas();
	vector<double> vegasR2 = P_RADmc2.GetVegas();
	vector<double> thetasR2 = P_RADmc2.GetThetas();
	vector<double> rhosR2 = P_RADmc2.GetRhos();
	for (int i = 0; i < a; ++i)
	{
		totalDeltaR2 += deltasR2[i];
		totalVegaR2 += vegasR2[i];
		totalThetaR2 += thetasR2[i];
		totalRhoR2 += rhosR2[i];
	}
	cout  << "Error ReverseAD MC2 total delta = " << abs(totalDelta - totalDeltaR2) << endl;
	cout  << "ReverseAD MC2 total vega = " << totalVegaR2 << endl;
	cout  << "ReverseAD MC2 total theta = " << totalThetaR2 << endl;
	cout  << "ReverseAD MC2 total rho = " << totalRhoR2 << endl;


	int bump = 1;
	vector<double>sgm(1, Sigma);
	vector<double>tt(1, t);
	vector<double>rr(1, r);
	double totalDeltafin = 0;	
	vector<double> deltasfin = P_DBL.FiniteDiff(S, bump);
	for (int i = 0; i < a; ++i)
	{
		totalDeltafin += deltasfin[i];
	}
	cout  << "Theory - FiniteDiff 0 total delta = " << abs(totalDelta - totalDeltafin) << endl;
	cout  << "Theory - FiniteDiff 0 total vega = " << abs(totalVega - P_DBL.FiniteDiff(sgm, bump)[0]) << endl;
	cout  << "Theory - FiniteDiff 0 total theta = " << abs(totalTheta - P_DBL.FiniteDiff(tt, bump)[0]) << endl;
	cout  << "Theory - FiniteDiff 0 total rho = " << abs(totalRho - P_DBL.FiniteDiff(tt, bump)[0]) << endl;
	
	double totalDeltafin1 = 0;
	vector<double> deltasfin1 = P_DBLmc1.FiniteDiff(S, bump);
	for (int i = 0; i < a; ++i)
	{
		totalDeltafin1 += deltasfin1[i];
	}
	cout << "Theory - FiniteDiff MC1 total delta = " << abs(totalDelta - totalDeltafin1) << endl;
	cout << "Theory - FiniteDiff MC1 total vega = " << abs(totalVega - P_DBLmc1.FiniteDiff(sgm, bump)[0]) << endl;
	cout << "Theory - FiniteDiff MC1 total theta = " << abs(totalTheta - P_DBLmc1.FiniteDiff(tt, bump)[0]) << endl;
	cout << "Theory - FiniteDiff MC1 total rho = " << abs(totalRho - P_DBLmc1.FiniteDiff(tt, bump)[0]) << endl;

	double totalDeltafin2 = 0;
	vector<double> deltasfin2 = P_DBLmc1.FiniteDiff(S, bump);
	for (int i = 0; i < a; ++i)
	{
		totalDeltafin2 += deltasfin2[i];
	}
	cout << "Theory - FiniteDiff MC2 total delta = " << abs(totalDelta - totalDeltafin2) << endl;
	cout << "Theory - FiniteDiff MC2 total vega = " << abs(totalVega - P_DBLmc2.FiniteDiff(sgm, bump)[0]) << endl;
	cout << "Theory - FiniteDiff MC2 total theta = " << abs(totalTheta - P_DBLmc2.FiniteDiff(tt, bump)[0]) << endl;
	cout << "Theory - FiniteDiff MC2 total rho = " << abs(totalRho - P_DBLmc2.FiniteDiff(tt, bump)[0]) << endl;

	time = clock();
	for (int i = 0; i < 1000*it; ++i)
	{
		P_DBL.FiniteDiff(S, bump);
		P_DBL.FiniteDiff(sgm, bump);
		P_DBL.FiniteDiff(tt, bump);
		P_DBL.FiniteDiff(rr, bump);
	}
	cout << "DBL FiniteDiff 0 time: " << double(clock() - time) / CLOCKS_PER_SEC / 1000 << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_DBLmc1.FiniteDiff(S, bump);
		P_DBLmc1.FiniteDiff(sgm, bump);
		P_DBLmc1.FiniteDiff(tt, bump);
		P_DBLmc1.FiniteDiff(rr, bump);
	}
	cout << "DBL FiniteDiff MC1 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;

	time = clock();
	for (int i = 0; i < it; ++i)
	{
		P_DBLmc2.FiniteDiff(S, bump);
		P_DBLmc2.FiniteDiff(sgm, bump);
		P_DBLmc2.FiniteDiff(tt, bump);
		P_DBLmc2.FiniteDiff(rr, bump);
	}
	cout << "DBL FiniteDiff MC2 time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;

	//cout << "------------0-------------" << endl;
	//cout << "DBL function value = " << P_DBL.GetPrice() << endl;
	//cout << "ReverseAD function value = " << P_RAD.GetPrice() << endl;
	//cout << "ForwardAD function value = " << P_FAD.GetPrice() << endl;
	//cout << "DBL function value - ReverseAD function value = " << abs(P_DBL.GetPrice()- P_RAD.GetPrice()) << endl;
	//cout << "DBL function value - ForwardAD function value = " << abs(P_DBL.GetPrice() - P_FAD.GetPrice()) << endl;
	//cout << "------------MC1-------------" << endl;
	//cout << "DBL function value = " << P_DBLmc1.GetPrice() << endl;
	//cout << "Theory - DBL function value = " << P_DBL.GetPrice() - P_DBLmc1.GetPrice() << endl;
	//cout << "ReverseAD function value = " << P_RADmc1.GetPrice() << endl;
	//cout << "ForwardAD function value = " << P_FADmc1.GetPrice() << endl;
	//cout << "Theory - ReverseAD function value = " << abs(P_DBL.GetPrice() - P_RADmc1.GetPrice()) << endl;
	//cout << "Theory - ForwardAD function value = " << abs(P_DBL.GetPrice() - P_FADmc1.GetPrice()) << endl;
	//cout << "------------MC2-------------" << endl;
	//cout << "DBL function value = " << P_DBLmc2.GetPrice() << endl;
	//cout << "Theory - DBL function value = " << P_DBL.GetPrice() - P_DBLmc2.GetPrice() << endl;
	//cout << "ReverseAD function value = " << P_RADmc2.GetPrice() << endl;
	//cout << "ForwardAD function value = " << P_FADmc2.GetPrice() << endl;
	//cout << "Theory - ReverseAD function value = " << abs(P_DBL.GetPrice() - P_RADmc2.GetPrice()) << endl;
	//cout << "Theory - ForwardAD function value = " << abs(P_DBL.GetPrice() - P_FADmc2.GetPrice()) << endl;

}

void SpotTest(const int it, const int n)
{
	double aS = 100;
	double bS = 500;
	vector<double> S = VarVector(aS, bS, n);
	vector<ForwardAD> SAD = VarVectorForwardAD(aS, bS, n); // will not work if we have other AD variables except S
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
	vector<double> S = VarVector(aS, bS, a); vector<ForwardAD>  SAD = VarVectorForwardAD(aS, bS, a);
	vector<double> sigma = VarVector(aSigma, bSigma, a); vector<ForwardAD> sigmReverseAD = VarVectorForwardAD(aSigma, bSigma, a);
	vector<double> 	r = VarVector(aR, bR, a); vector<ForwardAD>  rAD = VarVectorForwardAD(aR, bR, a);
	vector<double> t = VarVector(aT, aT, a); vector<ForwardAD> tAD = VarVectorForwardAD(aT, aT, a);
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
	vector<ForwardAD> SAD = VarVectorForwardAD(aS, bS, a); // will not work if we have other AD variables except S
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
	vector<ReverseAD> SAD1 = VarVectorReverseAD(aS, bS, a); 
	vector<ReverseAD> SAD2 = VarVectorReverseAD(aS, bS, a); 
	double K = 300;
	double Sigma(0.25);
	ReverseAD SigmaAD(Sigma);
	double r = 0.03;
	double T = 1;
	double t = 1. / 365;

	// compose 2 portfolios: one for double arguments and one for ForwardAD arguments
	vector<contract<double, double, double, double, double> > contractsDBL;
	vector<contract<ReverseAD, ReverseAD, double, double, ReverseAD> > contractsReverseAD;
	for (int i = 0; i < a; i++)
	{
		contract<double, double, double, double, double> cDBL(&S[i], &Sigma, &t, &r, T, K);
		contractsDBL.push_back(cDBL);
		contract<ReverseAD, ReverseAD, double, double, ReverseAD> cAD(&SAD1[i], &SigmaAD, &t, &r, T, K);
		contractsReverseAD.push_back(cAD);
	}

	portfolio<contract<double, double, double, double, double>, double> P_DBL(contractsDBL, "0");
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

	portfolio<contract<ReverseAD, ReverseAD, double, double, ReverseAD>, ReverseAD> P_ADmc1(contractsReverseAD, "MC1", n, m);
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

	portfolio<contract<ReverseAD, ReverseAD, double, double, ReverseAD>, ReverseAD> P_ADmc2(contractsReverseAD, "MC2", n, m);
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

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc1(contractsDBL, "MC1", n, m);
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

	portfolio<contract<double, double, double, double, double>, double> P_DBLmc2(contractsDBL, "MC2", n, m);
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


