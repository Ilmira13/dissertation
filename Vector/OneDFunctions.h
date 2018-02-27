#pragma once
#include <autodiff.h>

void FunctionsTests();

enum FunctionNames
{
	Function1,
	Function2,
	Function3
};

template <typename varType> class Functions
{
private:
	FunctionNames funcName;
	double bumpSize;
public:
	Functions(FunctionNames functionName = FunctionNames::Function1, double bmpSize = 0.0): funcName(functionName), bumpSize(bmpSize) {};

	double GetValue(double x)
	{
		if (typeid(varType) == typeid(autodiff))
		{
			vector<double> d(1, 1.0);
			autodiff newX(x, d);
			return Value(newX).value;
		}
		else
			return Value(x);
	}

	double GetDeriv(double x)
	{
		if (typeid(varType) == typeid(autodiff))
		{
			vector<double> d(1, 1.0);
			autodiff newX(x, d);
			return Value(newX).deriv[0];
		}
		else
		{
			double f1 = GetValue(x);
			double f2 = GetValue(x + bumpSize);

			return (f2 - f1) / bumpSize;
		}
	}

	template <typename varType> varType Value(varType x)
	{
		switch (funcName)
		{
		case FunctionNames::Function1:
			return Function1(x);
		case FunctionNames::Function2:
			return Function2(x);
		case FunctionNames::Function3:
			return Function3(x);
		default:
			throw exception("unknown function name: " + funcName);
		}
	}

	/*
	max(x, 0) + 100
	*/
	template <typename varType> varType Function1(varType x)
	{
		varType f = x;
		if (x <= 0.0)
			f = f - x;
		
		return f + 100;
	}

	/*
	1(x) + 100
	*/
	template <typename varType> varType Function2(varType x)
	{
		varType f = x - x + 1;
		if (x <= 0.0)
			f = f - 1;

		return f + 100;
	}

	/*
	ln(x)
	*/
	template <typename varType> varType Function3(varType x)
	{
		return log(x);
	}
};


