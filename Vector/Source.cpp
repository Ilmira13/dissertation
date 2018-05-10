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
	//SpotTest(1, 4);
	//FunctionsTests();
	//ReverseAD_tests();
	//ReverseADSpotTest(10, 40);
	//ReverseADSpotTest(100, 50); // 50 variables can show the difference
    //ReverseADMonteCarloTest(1, 1, 100, 1000); //const int it, const int a, const int n, const int m
	//dissertation();
	ReverseTest();
	system("pause");
	return 0;

}
