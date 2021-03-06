#include "pch.h"
#include "ForwardAD.h"

int i = 8;
double j = 7;
TEST(Constructor, DefaultConstructor) {
	ForwardAD c;
	vector<double> f;
	EXPECT_EQ(c.value, 0);
	EXPECT_EQ(c.deriv, f);
}

TEST(Constructor, DVDConstructor) {
	vector<double> f(i, j);
	ForwardAD c(j, f);
	EXPECT_EQ(c.value, j);
	EXPECT_EQ(size(c.deriv), i);
	for (int k = 0; k < i; k++)
	{
		EXPECT_EQ(c.deriv[k], j);
	}
}

TEST(Constructor, DIConstructor) {
	vector<double> f(i);
	ForwardAD c(j, f);
	EXPECT_EQ(c.value, j);
	EXPECT_EQ(size(c.deriv), i);
	for (int k = 0; k < i; k++)
	{
		EXPECT_EQ(c.deriv[k], 0);
	}
}

TEST(Constructor, DConstructor) {
	ForwardAD c;
	c.value = j;
	vector<double> f;
	EXPECT_EQ(c.value, j);
	EXPECT_EQ(c.deriv, f);
}

TEST(ClassOperator, EqualsOperatorForwardAD) {
	vector<double> f(i, j);
	ForwardAD b(double(i), f);
	ForwardAD c, a;
	c = a;	
	EXPECT_EQ(c.value, a.value);
	EXPECT_EQ(c.deriv, a.deriv);
	c = b;
	EXPECT_EQ(c.value, b.value);
	EXPECT_EQ(c.deriv, b.deriv);
}

TEST(ClassOperator, EqualsOperatorDouble) {
	ForwardAD c;
	c = j;
	EXPECT_EQ(c.value, j);
}

TEST(FriendVectorDoubleOperator, UnaryPlus) {
	vector<double> u(i, j), uu(i);
	uu = +u;
	for (int k = 0; k < i; k++)
	{
		EXPECT_EQ(uu[k], j);
	}
}

TEST(FriendVectorDoubleOperator, UnaryMinus) {
	vector<double> u(i, j), uu(i);
	uu = -u;
	for (int k = 0; k < i; k++)
	{
		EXPECT_EQ(uu[k], -j);
	}
}

TEST(FriendVectorDoubleOperator, Plus) {
	vector<double> u1(i, -j), u2(i+1, j), u(i), res;
	u.push_back(j);
	res = u1 + u2;
	EXPECT_EQ(size(u), size(res));
	for (int k = 0; k < i; k++)
	{
		EXPECT_EQ(u[k], res[k]);
	}
}