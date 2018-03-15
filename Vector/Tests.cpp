#include <Tests.h>

void FunctionsTests()
{
	std::cout.precision(8);

	double bumpSize = 0.0001;
	double x = 0.0;

	std::cout << "x = " << x << std::endl;
	std::cout << "Function1 = max(0.0, x) + 100" << std::endl;

	Functions<double> dblFunc1(FunctionNames::Function1, bumpSize);
	Functions<autodiff> adFunc1(FunctionNames::Function1);

	std::cout << "Function value = " << dblFunc1.GetValue(x) << std::endl;
	std::cout << "AD function value = " << adFunc1.GetValue(x) << std::endl;

	std::cout << "FD function derivative = " << dblFunc1.GetDeriv(x) << std::endl;
	std::cout << "AD function derivative = " << adFunc1.GetDeriv(x) << std::endl;

	std::cout << "Function2 = 1(x) + 100" << std::endl;

	Functions<double> dblFunc2(FunctionNames::Function2, bumpSize);
	Functions<autodiff> adFunc2(FunctionNames::Function2);

	std::cout << "Function value = " << dblFunc2.GetValue(x) << std::endl;
	std::cout << "AD function value = " << adFunc2.GetValue(x) << std::endl;

	std::cout << "FD function derivative = " << dblFunc2.GetDeriv(x) << std::endl;
	std::cout << "AD function derivative = " << adFunc2.GetDeriv(x) << std::endl;

	std::cout << "Function3 = ln(x)" << std::endl;

	Functions<double> dblFunc3(FunctionNames::Function3, bumpSize);
	Functions<autodiff> adFunc3(FunctionNames::Function3);

	std::cout << "Function value = " << dblFunc3.GetValue(x) << std::endl;
	std::cout << "AD function value = " << adFunc3.GetValue(x) << std::endl;

	std::cout << "FD function derivative = " << dblFunc3.GetDeriv(x) << std::endl;
	std::cout << "AD function derivative = " << adFunc3.GetDeriv(x) << std::endl;

	//system("pause");
}

void AAD_mult_test()
{
	// variables
	aad x1(1);
	aad x2(2);

	// calculations
	aad f = x1 * x2;
	//f.SetGradient(1.0);

	double x1g = x1.GetGradient();
	double x2g = x2.GetGradient();

	if (x1g == 2.0)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: multiplication test" << std::endl;
	if (x2g == 1.0)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: multiplication test" << std::endl;
}

void AAD_sum_test()
{
	// variables
	aad x1(1);
	aad x2(2);

	// calculations
	aad f = x1 + x2;
	//f.SetGradient(1.0);

	double x1g = x1.GetGradient();
	double x2g = x2.GetGradient();

	if (x1g == 1.0)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: summation test" << std::endl;
	if (x2g == 1.0)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: summation test" << std::endl;
}

void AAD_division_test()
{
	// variables
	aad x1(1);
	aad x2(2);

	// calculations
	aad f = x1 / x2;
	//f.SetGradient(1.0);

	double x1g = x1.GetGradient();
	double x2g = x2.GetGradient();

	if (x1g == 0.5)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: division test" << std::endl;
	if (x2g == -0.25)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: division test" << std::endl;
}

void AAD_vs_autodiff_test()
{
	// ----- aad ------------------
	// variables
	aad x1(1);
	aad x2(2);

	// calculations
	aad f0 = x1 / x2 + x1 * x1 * cos(x2 - x1);
	aad f1 = exp(f0);
	aad f2 = f1 + x1;
	//f2.SetGradient(1.0);

	double x1g = x1.GetGradient();
	double x2g = x2.GetGradient();

	// ----- autodiff --------------
	// variables
	autodiff x1a(1.0, 2);
	x1a.deriv[0] = 1.0;
	autodiff x2a(2.0, 2);
	x2a.deriv[1] = 1.0;

	// calculations
	autodiff f0a = x1a / x2a + x1a * x1a * cos(x2a - x1a);
	autodiff f1a = exp(f0a);
	autodiff f2a = f1a + x1a;

	// ----- comparison --------------
	if (x1g == f2a.deriv[0])
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: AAD vs autodiff test" << std::endl;
	if (x2g == f2a.deriv[1])
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: AAD vs autodiff test" << std::endl;
}

void AAD_several_runs_test()
{
	// variables
	aad x1(1);
	aad x2(2);

	// calculations
	aad f(0.0);
	int n = 2;
	for (int i = 0; i < n; ++i)
	{
		x1.Reset();
		x2.Reset();
		//f.Reset();
		f = x1 / x2;
	}
	//f.SetGradient(1.0);

	double x1g = x1.GetGradient();
	double x2g = x2.GetGradient();

	if (x1g == 0.5)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: several runs test" << std::endl;
	if (x2g == -0.25)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: several runs test" << std::endl;
}

void AAD_recursion_test()
{
	// variables
	aad x1(1);
	aad x2(2);

	// calculations
	aad f(0.0);
	int n = 10;
	for (int i = 0; i < n; ++i)
		f = f + x1 / x2;
	//f.SetGradient(1.0);

	double x1g = x1.GetGradient();
	double x2g = x2.GetGradient();

	if (x1g == 0.5 * n)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: recursion test" << std::endl;
	if (x2g == -0.25 * n)
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: recursion test" << std::endl;
}

void AAD_recursion_test2()
{
	// variables
	aad x1(1);
	double r = 0.1;

	// calculations
	int n = 10;
	for (int i = 1; i < n; ++i)
		x1 = x1 * r;
	//f.SetGradient(1.0);

	double x1g = x1.GetGradient();

	if (x1g == pow(r, n - 1))
		std::cout << "OK" << std::endl;
	else
		std::cout << "OK: recursion test #2 - attempt 1 failed" << std::endl; // this test should fail

	aad x2(1);
	vector<aad> f2(n, x2);
	for (int i = 1; i < n; ++i)
		f2[i] = f2[i - 1] * r;

	double x2g = f2[0].GetGradient();

	if (x2g == pow(r, n - 1))
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: recursion test #2 - attempt 2" << std::endl; // this test should pass
}

void AAD_many_variables_test()
{
	// variables
	int n = 4;
	aad x[4];
	x[0] = aad(1);
	x[1] = aad(2);
	x[2] = aad(3);
	x[3] = aad(4);

	// calculations
	aad f(0.0);
	for (int i = 0; i < n; ++i)
		f = f + x[i] * (n - i) + i;
	//f.SetGradient(1.0);

	for (int i = 0; i < n; ++i)
	{
		if (x[i].GetGradient() == (n - i))
			std::cout << "OK" << std::endl;
		else
			std::cout << "ERR: many variables test" << std::endl;
	}
}


void AAD_tests()
{
	AAD_mult_test();
	AAD_sum_test();
	AAD_division_test();
	AAD_vs_autodiff_test();

	AAD_several_runs_test();
	AAD_recursion_test();
	AAD_recursion_test2();
	AAD_many_variables_test();
}

