#pragma once
const double g = 9.8, pi = 3.14;
class simple_equations;

class simple_equations {
public:

	void count_d(double D, double delta_d)//считаем внутренний диаметр
	{
		diam_vnutr = (D - 2 * delta_d);
	}
	double get_d()//выводим внутренний диаметр
	{
		return diam_vnutr;
	}

	void count_epsilon(double delta)//считаем шероховатость
	{
		sherokh = delta/diam_vnutr;
	}
	double get_epsilon()//выводим шероховатость
	{
		return sherokh;
	}

	void count_v(double Q)//считаем скорость
	{
		speed = (4*Q)/(diam_vnutr * diam_vnutr * pi);
	}
	double get_v()//выводим скорость
	{
		return speed;
	}

	void count_Re(double nu)//считаем число Рейнольдса
	{
		Reynolds = (speed * diam_vnutr) / nu;
	}
	double get_Re()//выводим число Рейнольдса
	{
		return Reynolds;
	}

	void count_lambda()//считаем гидравлическое сопротивление
	{
		hydr_res = 0.11 * (pow((sherokh + (68.0 / Reynolds)), 0.25));
	}
	double get_lambda()//выводим гидравлическое сопротивление
	{
		return hydr_res;
	}

	void count_p_0(double L, double D, double delta_d, double delta,
		double z_0, double z_L, double rho, double nu, double p_n, double Q)//считаем давление на начале участка
	{
		count_d(D, delta_d);
		count_epsilon(delta);
		count_v(Q);
		count_Re(nu);
		count_lambda();
		p_c_0 = ((p_n / (rho * g)) + z_L - z_0 + (hydr_res * L * speed * speed) / (diam_vnutr * 2 * g))
			* (rho * g);
	}
	double get_p_0()//выводим давление на начале участка
	{
		return p_c_0;
	}

	void count_p_L(double L, double D, double delta_d, double delta,
		double z_0, double z_L, double rho, double nu, double p_n, double Q)//считаем давление на конце участка
	{
		count_d(D, delta_d);
		count_epsilon(delta);
		count_v(Q);
		count_Re(nu);
		count_lambda();
		p_c_L = ((p_n / (rho * g)) - z_L + z_0 - (hydr_res * L * speed * speed) / (diam_vnutr * 2 * g))
			* (rho * g);
	}
	double get_p_L()//выводим давление на конце участка
	{
		return p_c_L;
	}

	void count_Q(double L, double D, double delta_d, double delta,
		double z_0, double z_L, double rho, double nu, double p_0, double p_L)//считаем расход
	{
		count_d(D, delta_d);
		count_epsilon(delta);
		double lambda_approx;
		hydr_res = 0.02;
		double lamda_v_2 = (diam_vnutr * 2 * g * (((p_0 - p_L) / (rho * g)) + z_0 - z_L)) / L;
		do {
			lambda_approx = hydr_res;
			speed = sqrt(lamda_v_2 / lambda_approx);
			count_Re(nu);
			count_lambda();
		} while (abs(hydr_res - lambda_approx) > 0.0002);
		raskhod = pi * diam_vnutr * diam_vnutr * speed / 4;
	}
	double get_Q()//выводим расход
	{
		return raskhod;
	}

	//void subsequent_approaches()
	//{

	//}
private:
	double diam_vnutr, sherokh, speed, Reynolds, hydr_res, p_c_0, p_c_L, raskhod;
};

TEST(MOC_Solver, Task_1)
{
	//начальные условия
	double L = 80e3, D = 0.72, delta_d = 0.01,
		delta = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_L = 0.6e6, Q = 0.972;
	simple_equations simple;
	simple.count_p_0(L, D, delta_d, delta, z_0, z_L, rho, nu, p_L, Q);
	double p_0 = simple.get_p_0();
}

TEST(MOC_Solver, Task_2)
{
	//начальные условия
	double L = 80e3, D = 0.72, delta_d = 0.01,
		delta = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 5e6, p_L = 0.8e6;
	simple_equations simple;
	simple.count_Q(L, D, delta_d, delta, z_0, z_L, rho, nu, p_0, p_L);
	double Q = simple.get_Q();
}