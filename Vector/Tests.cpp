#include <Tests.h>

void FunctionsTests()
{
	std::cout.precision(8);

	double bumpSize = 0.0001;
	double x = 0.0;

	std::cout << "x = " << x << std::endl;
	std::cout << "Function1 = max(0.0, x) + 100" << std::endl;

	Functions<double> dblFunc1(FunctionNames::Function1, bumpSize);
	Functions<ForwardAD> adFunc1(FunctionNames::Function1);

	std::cout << "Function value = " << dblFunc1.GetValue(x) << std::endl;
	std::cout << "AD function value = " << adFunc1.GetValue(x) << std::endl;

	std::cout << "FD function derivative = " << dblFunc1.GetDeriv(x) << std::endl;
	std::cout << "AD function derivative = " << adFunc1.GetDeriv(x) << std::endl;

	std::cout << "Function2 = 1(x) + 100" << std::endl;

	Functions<double> dblFunc2(FunctionNames::Function2, bumpSize);
	Functions<ForwardAD> adFunc2(FunctionNames::Function2);

	std::cout << "Function value = " << dblFunc2.GetValue(x) << std::endl;
	std::cout << "AD function value = " << adFunc2.GetValue(x) << std::endl;

	std::cout << "FD function derivative = " << dblFunc2.GetDeriv(x) << std::endl;
	std::cout << "AD function derivative = " << adFunc2.GetDeriv(x) << std::endl;

	std::cout << "Function3 = ln(x)" << std::endl;

	Functions<double> dblFunc3(FunctionNames::Function3, bumpSize);
	Functions<ForwardAD> adFunc3(FunctionNames::Function3);

	std::cout << "Function value = " << dblFunc3.GetValue(x) << std::endl;
	std::cout << "AD function value = " << adFunc3.GetValue(x) << std::endl;

	std::cout << "FD function derivative = " << dblFunc3.GetDeriv(x) << std::endl;
	std::cout << "AD function derivative = " << adFunc3.GetDeriv(x) << std::endl;

	//system("pause");
}

void ReverseAD_mult_test()
{
	// variables
	ReverseAD x1(1);
	ReverseAD x2(2);

	// calculations
	ReverseAD f = x1 * x2;
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

void ReverseAD_sum_test()
{
	// variables
	ReverseAD x1(1);
	ReverseAD x2(2);

	// calculations
	ReverseAD f = x1 + x2;
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

void ReverseAD_division_test()
{
	// variables
	ReverseAD x1(1);
	ReverseAD x2(2);

	// calculations
	ReverseAD f = x1 / x2;
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

void ReverseAD_vs_ForwardAD_test()
{
	// ----- ReverseAD ------------------
	// variables
	ReverseAD x1(1);
	ReverseAD x2(2);

	// calculations
	ReverseAD f0 = x1 / x2 + x1 * x1 * cos(x2 - x1);
	ReverseAD f1 = exp(f0);
	ReverseAD f2 = f1 + x1;
	//f2.SetGradient(1.0);

	double x1g = x1.GetGradient();
	double x2g = x2.GetGradient();

	// ----- ForwardAD --------------
	// variables
	ForwardAD x1a(1.0, 2);
	x1a.deriv[0] = 1.0;
	ForwardAD x2a(2.0, 2);
	x2a.deriv[1] = 1.0;

	// calculations
	ForwardAD f0a = x1a / x2a + x1a * x1a * cos(x2a - x1a);
	ForwardAD f1a = exp(f0a);
	ForwardAD f2a = f1a + x1a;

	// ----- comparison --------------
	if (x1g == f2a.deriv[0])
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: ReverseAD vs ForwardAD test" << std::endl;
	if (x2g == f2a.deriv[1])
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: ReverseAD vs ForwardAD test" << std::endl;
}

void ReverseAD_several_runs_test()
{
	// variables
	ReverseAD x1(1);
	ReverseAD x2(2);

	// calculations
	ReverseAD f(0.0);
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

void ReverseAD_recursion_test()
{
	// variables
	ReverseAD x1(1);
	ReverseAD x2(2);

	// calculations
	ReverseAD f(0.0);
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

void ReverseAD_recursion_test2()
{
	// variables
	ReverseAD x1(1);
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

	ReverseAD x2(1);
	vector<ReverseAD> f2(n, x2);
	for (int i = 1; i < n; ++i)
		f2[i] = f2[i - 1] * r;

	double x2g = f2[0].GetGradient();

	if (x2g == pow(r, n - 1))
		std::cout << "OK" << std::endl;
	else
		std::cout << "ERR: recursion test #2 - attempt 2" << std::endl; // this test should pass
}

void ReverseAD_many_variables_test()
{
	// variables
	int n = 4;
	ReverseAD x[4];
	x[0] = ReverseAD(1);
	x[1] = ReverseAD(2);
	x[2] = ReverseAD(3);
	x[3] = ReverseAD(4);

	// calculations
	ReverseAD f(0.0);
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


void ReverseAD_tests()
{
	ReverseAD_mult_test();
	ReverseAD_sum_test();
	ReverseAD_division_test();
	ReverseAD_vs_ForwardAD_test();

	ReverseAD_several_runs_test();
	ReverseAD_recursion_test();
	ReverseAD_recursion_test2();
	ReverseAD_many_variables_test();
}

