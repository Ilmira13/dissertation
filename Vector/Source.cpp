#include <iostream>
#include <RND.h>
#include <BlackScholes.h>
#include <Portfolio.h>

using namespace std;

int main() {
	cout.precision(8);

	//MonteCarloTest();
	cout << "AllVarTest---------------------------------------------------------" << endl;
	AllVarTest(100, 4);
	cout << "SpotTest---------------------------------------------------------" << endl;
	SpotTest(100, 4);
	//AADtest();
	system("pause");
	return 0;

}
