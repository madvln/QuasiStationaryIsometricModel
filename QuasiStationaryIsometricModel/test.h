#pragma once
const double g = 9.81, pi = 3.14;
class simple_equations;
using namespace std;

class simple_equations : public fixed_system_t<1>
{
	double& L;
	double& D;
	double& delta_d;
	double& delta;
	double& z_0;
	double& z_L;
	double& rho;
	double& nu;
	double& p_0;
	double& p_L;
	double& Q;
	double& h;
	using fixed_system_t<1>::var_type;
public:
	simple_equations(double& L, double& D, double& delta_d, double& delta,
		double& z_0, double& z_L, double& rho, double& nu, double& p_0, double& p_L, double& Q, double& h) :
		L{ L }, D{ D }, delta_d{ delta_d }, delta{ delta }, z_0{ z_0 }, z_L{ z_L }, rho{ rho }, nu{ nu }, p_0{ p_0 }, p_L{ p_L }, Q{ Q }, h{ h }
	{
		n = static_cast<int>(L / h + 0.5) + 1;
	}
	// Задание функции невязок
	var_type residuals(const var_type& x) {
		speed = x;
		count_d();
		count_epsilon();
		count_Re();
		count_lambda();
		return
		{
			hydr_res * (L * speed * speed / (diam_vnutr * 2 * g)) + p_c.back() / (rho * g) + z_L - p_c[0] / (rho * g) + z_0
		};
	}

	void count_d()//считаем внутренний диаметр
	{
		diam_vnutr = D - 2 * delta_d;
	}
	double get_d()//выводим внутренний диаметр
	{
		return diam_vnutr;
	}

	void count_epsilon()//считаем шероховатость
	{
		sherokh = delta / diam_vnutr;
	}
	double get_epsilon()//выводим шероховатость
	{
		return sherokh;
	}

	void count_v()//считаем скорость
	{
		speed = (4 * Q) / (diam_vnutr * diam_vnutr * pi);
	}
	double get_v()//выводим скорость
	{
		return speed;
	}

	void count_Re()//считаем число Рейнольдса
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

	void count_p_0()//считаем давление на начале участка
	{
		count_d();
		count_epsilon();
		count_v();
		count_Re();
		count_lambda();
		p_c_0 = ((p_L / (rho * g)) + z_L - z_0 + (hydr_res * L * speed * speed) / (diam_vnutr * 2 * g))
			* (rho * g);
	}
	double get_p_0()//выводим давление на начале участка
	{
		return p_c_0;
	}

	void count_p_L()//считаем давление на конце участка
	{
		count_d();
		count_epsilon();
		count_v();
		count_Re();
		count_lambda();
		p_c_L = ((p_0 / (rho * g)) - z_L + z_0 - (hydr_res * L * speed * speed) / (diam_vnutr * 2 * g))
			* (rho * g);
	}
	double get_p_L()//выводим давление на конце участка
	{
		return p_c_L;
	}

	void count_Q()//считаем расход
	{
		count_d();
		count_epsilon();
		double lambda_approx;
		hydr_res = 0.02;
		double lamda_v_2 = (diam_vnutr * 2 * g * (((p_0 - p_L) / (rho * g)) + z_0 - z_L)) / L;
		do {
			lambda_approx = hydr_res;
			speed = sqrt(lamda_v_2 / lambda_approx);
			count_Re();
			count_lambda();
		} while (abs(hydr_res - lambda_approx) > 0.0002);
		raskhod = pi * diam_vnutr * diam_vnutr * speed / 4;
	}
	double get_Q()//выводим расход
	{
		return raskhod;
	}

	void count_tau()
	{
		count_d();
		count_epsilon();
		count_v();
		count_Re();
		count_lambda();
		tau = (hydr_res / 8) * rho * speed * speed;
	}

	void eiler_from_start()
	{
		double delta_z = (z_L - z_0) / (n - 1);
		p_c = vector<double>(n);
		count_tau();
		p_c[0] = p_0;
		for (int i = 1; i < n; i++)
			p_c[i] = p_c[i - 1] + h * ((-4 / diam_vnutr) * tau - rho * g * (delta_z / h));
	}

	

	void eiler_from_end()
	{
		double delta_z = (z_L - z_0) / (n - 1);
		p_c = vector<double>(n);
		count_tau();
		p_c[n - 1] = p_L;
		for (int i = n - 2; i >= 0; i--)
			p_c[i] = p_c[i + 1] - h * ((-4 / diam_vnutr) * tau - rho * g * (delta_z / h));
	}

	void newtone_from_end()
	{
		double delta_z = (z_L - z_0) / (n - 1);
		p_c = vector<double>(n);
		count_tau();
		p_c[n - 1] = p_L;
		for (int i = n - 2; i >= 0; i--)
			p_c[i] = p_c[i + 1] - ((-4 / diam_vnutr) * tau - rho * g * (delta_z / h)) * ((h * (i + 1)) - (h * i));
	}
	vector<double> get_p_profile()
	{
		return p_c;
	}

	int get_point_count() const {
		return n;
	}

private:
	int n;
	double diam_vnutr, sherokh, speed, Reynolds, hydr_res, p_c_0, p_c_L, raskhod, tau;
	vector<double> p_c;
};

TEST(MOC_Solver, Task_1)
{
	//начальные условия
	double L = 80e3, D = 0.72, delta_d = 0.01,
		delta = 15e-6, z_0 = 100, z_L = 50,
		rho = 870, nu = 15e-6, p_L = 0.6e6, Q = 0.972, h = 1e3;
	double p_0 = 0.0;
	simple_equations simple(L, D, delta_d, delta, z_0, z_L, rho, nu, p_0, p_L, Q, h);
	simple.count_p_0();
	p_0 = simple.get_p_0();
}

TEST(MOC_Solver, Task_2)
{
	//начальные условия
	double L = 80e3, D = 0.72, delta_d = 0.01,
		delta = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 5e6, p_L = 0.8e6, h = 1e3;
	double Q = 0.0;
	simple_equations simple(L, D, delta_d, delta, z_0, z_L, rho, nu, p_0, p_L, Q, h);
	simple.count_Q();
	Q = simple.get_Q()*3600.0;
}

// Функция, описывающая правую часть дифференциального уравнения

TEST(MOC_Solver, Task_3)
{
	double p_0 = 0.0;
	double p_L = 0.6e6;
	double L = 80e3; //Длина участка трубы
	double h = 1e3; //Шаг
	double D = 0.72; //внешний диаметр трубы
	double delta_d = 0.01; //толщина стенки
	double delta = 15e-6;
	double Q = 3500.0 / 3600.0;
	double rho = 870;
	double nu = 15e-6;
	double z_0 = 100;
	double z_L = 50;
	simple_equations simple(L, D, delta_d, delta, z_0, z_L, rho, nu, p_0, p_L, Q, h); //объявляем переменную класса для расчетов
	int n = simple.get_point_count(); //Количество шагов
	simple.eiler_from_end();
	vector<double> p_profile = simple.get_p_profile();
	simple.newtone_from_end();
	vector<double> p_profile_2 = simple.get_p_profile();
}


TEST(MOC_Solver, Task_4)
{
	double p_0 = 5.65e6;
	double p_L = 0.6e6;
	double L = 80e3; //Длина участка трубы
	double h = 1e3; //Шаг
	double D = 0.72; //внешний диаметр трубы
	double delta_d = 0.01; //толщина стенки
	double delta = 15e-6;
	double Q = 3500.0 / 3600.0;
	double rho = 870;
	double nu = 15e-6;
	double z_0 = 100;
	double z_L = 50;
	simple_equations simple(L, D, delta_d, delta, z_0, z_L, rho, nu, p_0, p_L, Q, h); //объявляем переменную класса для расчетов
	int n = simple.get_point_count(); //Количество шагов
	simple.newtone_from_end();
	vector<double> p_profile_2 = simple.get_p_profile();
	// Задание настроек решателя по умолчанию
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	fixed_newton_raphson<1>::solve_dense(simple, { 10 }, parameters, &result);
	double d = simple.get_d();
	cout << result.argument * pi * pow(d, 2) / 4 * 3600 << endl;
}