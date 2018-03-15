#include <iostream>
#include <Tests.h>
#include <RND.h>
#include <BlackScholes.h>
#include <Portfolio.h>

using namespace std;

int main() {
	cout.precision(8);

	//MonteCarloTest();
	//cout << "AllVarTest---------------------------------------------------------" << endl;
	//AllVarTest(100, 4);
	//cout << "SpotTest---------------------------------------------------------" << endl;
	SpotTest(100, 4);
	//FunctionsTests();
	//AAD_tests();
	//AADSpotTest(10, 40);
	//AADSpotTest(1, 1);
	//AADSpotTest(100, 50); // 50 variables can show the difference
	//AADMonteCarloTest(1, 1, 100, 1000); //const int it, const int a, const int n, const int m
	system("pause");
	return 0;

}
