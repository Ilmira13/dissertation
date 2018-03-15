#include "pch.h"
#include "autodiff.h"

TEST(TestCaseName, TestName) {
	vector<double> u;
    autodiff c(1, u);
  EXPECT_EQ(c.value, 0);

}