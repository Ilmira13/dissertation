#pragma once
#include <iostream>
#include <BlackScholes.h>
#include <Portfolio.h>
#include <ReverseAD.h>
#include <ctime>

using namespace std;

// contract <S, sigma, t, r, T, K>
//void SpotTest(const int it = 100, const int n = 1);
//void AllVarTest(const int it = 10, const int n = 4);
//void MonteCarloTest(const int it = 1, const int a = 1, const int n = 10000, const int m = 10000);
//void ReverseADSpotTest(const int it, const int n); // works for it = 1 only
//void ReverseADMonteCarloTest(const int it = 1, const int a = 1, const int n = 10, const int m = 1000); // does not work
void dissertation(const int it = 1, const int a = 10, const int n = 5, const int m = 100);
