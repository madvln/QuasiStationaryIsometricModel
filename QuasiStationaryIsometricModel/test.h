#pragma once
const double g = 9.81, pi = 3.14;
using namespace std;

struct my_pipe_parameters
{
	/// @brief Длина трубы, (м)
	double length = 80e3;
	/// @brief  Внешний диаметр трубы, (м)
	double D = 0.72;
	/// @brief Толщина стенки трубы, (м)
	double delta_d = 0.01;
	/// @brief Абсолютная шероховатость, (м)
	double abs_roughness = 15e-6;
	/// @brief Внутренний диаметр трубы, (м)
	double d = D - 2 * delta_d;
	/// @brief Шероховатость
	double roughness = abs_roughness / d;

	double z_0 = 100;

	double z_L = 50;

	my_pipe_parameters(double length, double D, double delta_d, double abs_roughness, double z_0, double z_L):
		length{ length }, D{ D }, delta_d { delta_d }, abs_roughness { abs_roughness }, z_0 { z_0 }, z_L { z_L }
	{
	}
};

struct my_task_parameters
{
	my_pipe_parameters& pipe;
	/// @brief Плотность жидкости, (кг/м3)
	double rho = 870;
	/// @brief Кинематическая вязкость, (м^2/с)
	double nu = 15e-6;
	/// @brief Давление в начале участка, (Па)
	double p_0 = 5e6;
	/// @brief Давление в конце участка, (Па)
	double p_L = 0.6e6;
	/// @brief Расход жидкости, (м3/с)
	double Q = 0.972;
	/// @brief Шаг сетки, (м)
	double h = 1e3;
	/// @brief Количество шагов
	double n = static_cast<int>(pipe.length / h + 0.5) + 1;
	/// @brief Профиль давлений
	vector<double> p_profile = vector<double>(n);
	my_task_parameters(my_pipe_parameters& pipe, double rho, double nu, double p_0, double p_L, double Q, double h) :
		pipe{ pipe }, rho{ rho }, nu{ nu }, p_0{ p_0 }, p_L{ p_L }, Q{ Q }, h{ h }
	{		
	}
};

/// @brief Функция расчета скорости из расхода
/// @param Q Расход
/// @param d Внутренний диаметр
/// @return Скорость
double count_speed(double Q, double d)
{
	double speed = (4 * Q) / (pow(d, 2) * pi);
	return speed;
}

/// @brief Функция расчета числа Рейнольдса
/// @param speed Скорость, (м/с)
/// @param d Внутренний диаметр, (м)
/// @param nu Кинетическая вязкость, (м^2/с)
/// @return Число Рейнольдса
double count_Re(double speed, double d, double nu)
{
	double Re = (speed * d) / nu;
	return Re;
}

/// @brief Функция расчета лямбды с помощью pde_solvers.hydraulic_resistance_isaev
/// @param Re Число Рейнольдса
/// @param roughness Шероховатость
/// @return Лямбда
double count_lambda(double Re, double roughness)
{
	double lambda = hydraulic_resistance_isaev(Re, roughness);
	return lambda;
}

double count_tau(double lambda, double rho, double speed)
{
	double tau = (lambda / 8) * rho * pow(speed, 2);
	return tau;
}

class classic_solver
{
	my_pipe_parameters& pipe;
	my_task_parameters& task;
public:
	classic_solver(my_pipe_parameters& pipe, my_task_parameters& task) :
		pipe{ pipe }, task{ task }
	{
	}
	double QP_task()
	{
		/// @brief Скорость, (м/с)
		double speed = count_speed(task.Q, pipe.d);
		/// @brief Число Рейнольдса
		double Re = count_Re(speed, pipe.d, task.nu);
		/// @brief Лямбда
		double lambda = count_lambda(Re, pipe.roughness);
		double p_0 = task.p_L + (pipe.z_0 - pipe.z_L) * task.rho * g +
			(lambda * pipe.length * pow(speed, 2) * task.rho) / (2 * task.pipe.d);
		return p_0;
	}
	double PQ_task()
	{
		/// @brief Скорость, (м/с)
		double speed = count_speed(task.Q, pipe.d);
		/// @brief Число Рейнольдса
		double Re = count_Re(speed, pipe.d, task.nu);
		/// @brief Лямбда
		double lambda = count_lambda(Re, pipe.roughness);
		double p_L = task.p_0 - (pipe.z_0 - pipe.z_L) * task.rho * g -
			(lambda * pipe.length * pow(speed, 2) * task.rho) / (2 * task.pipe.d);
		return p_L;
	}
	double PP_task()
	{
		double lambda_prev, speed, Re;
		double lambda = 0.02;
		double lambda_v2 = (pipe.d * 2 * g * (((task.p_0 - task.p_L) / (task.rho * g)) + pipe.z_0 - pipe.z_L)) / pipe.length;
		do {
			lambda_prev = lambda;
			speed = sqrt(lambda_v2 / lambda_prev);
			Re = count_Re(speed, pipe.d, task.nu);
			lambda = count_lambda(Re, pipe.roughness);
		} while (abs(lambda - lambda_prev) > 0.0002);
		double Q = pi * pow(pipe.d, 2) * speed / 4;
		return Q;
	}
};

class euler_solver : public fixed_system_t<1>
{
	my_pipe_parameters& pipe;
	my_task_parameters& task;
	using fixed_system_t<1>::var_type;
public:
	euler_solver(my_pipe_parameters& pipe, my_task_parameters& task) :
		pipe{ pipe }, task{ task }
	{
	}

	var_type residuals(const var_type& x) 
	{
		double speed = x;
		double Re = count_Re(speed, pipe.d, task.nu);
		double lambda = count_lambda(Re, pipe.roughness);
		return (lambda * (pipe.length * pow(speed, 2) / (pipe.d * 2 * g)) + task.p_L / (task.rho * g) + pipe.z_L - task.p_0 / (task.rho * g) - pipe.z_0);
	}

	vector<double> euler_from_start()
	{
		double delta_z = (pipe.z_0 - pipe.z_L) / (task.n);
		double speed = count_speed(task.Q, pipe.d);
		double Re = count_Re(speed, pipe.d, task.nu);
		double lambda = count_lambda(Re, pipe.roughness);
		double tau = count_tau(lambda, task.rho, speed);
		task.p_profile[0] = task.p_0;
		for (int i = 1; i < task.n; i++)
			task.p_profile[i] = task.p_profile[i - 1] + task.h * ((-4 / pipe.d) * tau - task.rho * g * (delta_z / task.h));
		task.p_L = task.p_profile.back();
		return task.p_profile;
	}

	vector<double> euler_from_end()
	{
		double delta_z = (pipe.z_0 - pipe.z_L) / (task.n - 1);
		double speed = count_speed(task.Q, pipe.d);
		double Re = count_Re(speed, pipe.d, task.nu);
		double lambda = count_lambda(Re, pipe.roughness);
		double tau = count_tau(lambda, task.rho, speed);
		task.p_profile[task.n - 1] = task.p_L;
		for (int i = task.n - 2; i >= 0; i--)
			task.p_profile[i] = task.p_profile[i + 1] - task.h * ((-4 / pipe.d) * tau - task.rho * g * (delta_z / task.h));
		task.p_0 = task.p_profile.back();
		return task.p_profile;
	}
};

TEST(MOC_Solver, Task_1)
{
	simple_pipe_properties simple_pipe;
	/// @brief 
	simple_pipe.diameter = 0.72;
	/// @brief 
	simple_pipe.length = 80e3;
	/// @brief 
	simple_pipe.dx = 1e3;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 0.0, p_L = 0.6e6, Q = 0.972;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, simple_pipe.dx };
	classic_solver simple(pipe,task);
	p_0 = simple.QP_task();
}

TEST(MOC_Solver, Task_2)
{
	simple_pipe_properties simple_pipe;
	/// @brief 
	simple_pipe.diameter = 0.72;
	/// @brief 
	simple_pipe.length = 80e3;
	/// @brief 
	simple_pipe.dx = 1e3;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 5e6, p_L = 0.8e6, Q = 0.0;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, simple_pipe.dx };
	classic_solver simple(pipe, task);
	Q = simple.PP_task() * 3600;
}

TEST(MOC_Solver, Task_3)
{
	simple_pipe_properties simple_pipe;
	/// @brief 
	simple_pipe.diameter = 0.72;
	/// @brief 
	simple_pipe.length = 80e3;
	/// @brief 
	simple_pipe.dx = 1e3;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 0.0, p_L = 0.6e6, Q = 0.972;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, simple_pipe.dx };
	euler_solver e_solver(pipe, task);
	vector<double> p_profile = e_solver.euler_from_end();
}

TEST(MOC_Solver, Task_4)
{
	simple_pipe_properties simple_pipe;
	/// @brief 
	simple_pipe.diameter = 0.72;
	/// @brief 
	simple_pipe.length = 80e3;
	/// @brief 
	simple_pipe.dx = 1e3;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 5e6, p_L = 0.8e6, Q = 0.0;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, simple_pipe.dx };
	euler_solver e_solver(pipe, task);
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	fixed_newton_raphson<1>::solve_dense(e_solver, { 20 }, parameters, &result);
	cout << result.argument * pi * pow(pipe.d, 2) / 4 * 3600 << endl;
}