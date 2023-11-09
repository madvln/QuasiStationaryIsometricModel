#pragma once
const double g = 9.8, pi = 3.14;
class simple_equations;

struct simple_struct {
	double L; double D; double delta_d; double delta;
	double z_0; double z_L; double rho; double nu; double p_n; double Q; double p_0; double p_L;
};

class simple_equations {
public:

	void count_d(double D, double delta_d)//������� ���������� �������
	{
		diam_vnutr = (D - 2 * delta_d);
	}
	double get_d()//������� ���������� �������
	{
		return diam_vnutr;
	}

	void count_epsilon(double delta)//������� �������������
	{
		sherokh = delta / diam_vnutr;
	}
	double get_epsilon()//������� �������������
	{
		return sherokh;
	}

	void count_v(double Q)//������� ��������
	{
		speed = (4 * Q) / (diam_vnutr * diam_vnutr * pi);
	}
	double get_v()//������� ��������
	{
		return speed;
	}

	void count_Re(double nu)//������� ����� ����������
	{
		Reynolds = (speed * diam_vnutr) / nu;
	}
	double get_Re()//������� ����� ����������
	{
		return Reynolds;
	}

	void count_lambda()//������� �������������� �������������
	{
		hydr_res = 0.11 * (pow((sherokh + (68.0 / Reynolds)), 0.25));
	}
	double get_lambda()//������� �������������� �������������
	{
		return hydr_res;
	}

	void count_p_0(double L, double D, double delta_d, double delta,
		double z_0, double z_L, double rho, double nu, double p_n, double Q)//считаем давление на начале участка
	{
		
		count_d(pimple_.D, pimple_.delta_d);
		count_epsilon(pimple_.delta);
		count_v(pimple_.Q);
		count_Re(pimple_.nu);
		count_lambda();
		p_c_0 = ((p_n / (rho * g)) + z_L - z_0 + (hydr_res * L * speed * speed) / (diam_vnutr * 2 * g))
			* (rho * g);
	}
	double get_p_0()
	{
		return p_c_0;
	}

	void count_p_L(double L, double z_0, double z_L, double rho, double p_n)
	{
		count_d(D, delta_d);
		count_epsilon(delta);
		count_v(Q);
		count_Re(nu);
		count_lambda();
		p_c_L = ((p_n / (rho * g)) - z_L + z_0 - (hydr_res * L * speed * speed) / (diam_vnutr * 2 * g))
			* (rho * g);
	}
	double get_p_L()//������� �������� �� ����� �������
	{
		return p_c_L;
	}

	void count_Q(double L, double D, double delta_d, double delta,
		double z_0, double z_L, double rho, double nu, double p_0, double p_L)//������� ������
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
	double get_Q()//������� ������
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
	//��������� �������
	double L = 80e3, D = 0.72, delta_d = 0.01,
		delta = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_L = 0.6e6, Q = 0.972, p_0;
	simple_equations simple;
	simple_struct pimple{ L,D,delta_d, delta, z_0, z_L, rho, nu, Q, p_0, p_L};
	simple.count_p_0(pimple);
	p_0 = simple.get_p_0();
}

TEST(MOC_Solver, Task_2)
{
	//��������� �������
	double L = 80e3, D = 0.72, delta_d = 0.01,
		delta = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, Q, p_0 = 5e6, p_L = 0.8e6;
	simple_equations simple;
	simple_struct pimple{L,D,delta_d, delta, z_0, z_L, rho, nu, Q, p_0, p_L};
	pimple.D;
	simple.count_Q(L, D, delta_d, delta, z_0, z_L, rho, nu, p_0, p_L);
	Q = simple.get_Q();
}