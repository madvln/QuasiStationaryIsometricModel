#pragma once

class simple_equations;

class simple_equations {
public:
	simple_equations(double& D, double& delta_d)
		:D_(D), delta_d_(delta_d)
	{
		count_d(D, delta_d);
	}

	void count_d(double D, double delta_d)
	{
		diam_vnutr = D - 2 * delta_d;
	}
	double get_d()
	{
		return diam_vnutr;
	}


private:
	double& D_, delta_d_;
	double diam_vnutr;
};

TEST(MOC_Solver, Task_1)
{
	const double g = 9.8;
	//начальные условия
	double D = 720, delta_d = 10,
		delta = 0.015, z_0 = 50, z_L = 100,
		rho = 870, nu = 15, p_L = 0.6, Q = 3500;
	//double d = D - 2 * delta_d;
	simple_equations simple(D, delta_d);
	double d = simple.get_d();

}
