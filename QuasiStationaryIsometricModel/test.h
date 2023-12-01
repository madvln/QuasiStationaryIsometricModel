#pragma once
const double g = 9.81, pi = 3.14;
using namespace std;

struct my_pipe_parameters
{
	/// @brief Длина трубы, (м)
	double length = 80e3;
	/// @brief  Внешний диаметр трубы, (м)
	double external_diameter = 0.72;
	/// @brief Толщина стенки трубы, (м)
	double delta_d = 0.01;
	/// @brief Абсолютная шероховатость, (м)
	double abs_roughness = 15e-6;
	/// @brief Внутренний диаметр трубы, (м)
	double internal_diameter = external_diameter - 2 * delta_d;
	/// @brief Шероховатость
	double roughness = abs_roughness / internal_diameter;
	/// @brief Начальная высотная отметка, (м)
	double z_0 = 100;
	/// @brief Конечная высотная отметка, (м)
	double z_L = 50;
	/// @brief Шаг сетки, (м)
	double h = 1e3;
	/// @brief Количество шагов
	size_t n;
	my_pipe_parameters(double length, double external_diameter, double delta_d, double abs_roughness, double z_0, double z_L, double h):
		length{ length }, external_diameter{ external_diameter }, delta_d { delta_d }, abs_roughness { abs_roughness }, z_0 { z_0 }, z_L { z_L }, h{ h }
	{
		n = static_cast<int>(length / h + 0.5) + 1;
	}
};

struct my_task_parameters
{
	my_pipe_parameters& pipe;
	/// @brief Плотность жидкости, (кг/м^3)
	double rho = 870;
	/// @brief Кинематическая вязкость, (м^2/с)
	double nu = 15e-6;
	/// @brief Давление в начале участка, (Па)
	double p_0 = 5e6;
	/// @brief Давление в конце участка, (Па)
	double p_L = 0.6e6;
	/// @brief Расход жидкости, (м^3/с)
	double Q = 0.972;
	my_task_parameters(my_pipe_parameters& pipe, double rho, double nu, double p_0, double p_L, double Q) :
		pipe{ pipe }, rho{ rho }, nu{ nu }, p_0{ p_0 }, p_L{ p_L }, Q{ Q }
	{	
	}
};

/// @brief Функция расчета скорости из расхода
/// @param Q Расход, (м^3/с)
/// @param internal_diameter Внутренний диаметр, (м)
/// @return Скорость, (м/с)
double calc_speed(double Q, double internal_diameter)
{
	double speed = (4 * Q) / (pow(internal_diameter, 2) * pi);
	return speed;
}

/// @brief Функция расчета числа Рейнольдса
/// @param speed Скорость, (м/с)
/// @param internal_diameter Внутренний диаметр, (м)
/// @param nu Кинетическая вязкость, (м^2/с)
/// @return Число Рейнольдса
double calc_Re(double speed, double internal_diameter, double nu)
{
	double Re = (speed * internal_diameter) / nu;
	return Re;
}

/// @brief Функция расчета касательного напряжения трения
/// @param hydraulic_resistance гидравлическое сопротивление
/// @param rho Плотность, (кг/см^2)
/// @param speed Скорость, (м/с)
/// @return касательное напряжение трения
double calc_tau(double hydraulic_resistance, double rho, double speed)
{
	double tau = (hydraulic_resistance / 8) * rho * pow(speed, 2);
	return tau;
}

/// @brief класс, решающий задачи QP, PQ, PP методом простых итераций
class simple_iterations_solver
{
	const my_pipe_parameters& pipe;
	const my_task_parameters& task;
public:
	simple_iterations_solver(my_pipe_parameters& pipe, my_task_parameters& task) :
		pipe(pipe), task(task)
	{
	}
	/// @brief Метод нахождения входного давления
	/// @return Pвх
	double QP_task()
	{
		double speed = calc_speed(task.Q, pipe.internal_diameter);
		double Re = calc_Re(speed, pipe.internal_diameter, task.nu);
		double hydraulic_resistance = hydraulic_resistance_isaev(Re, pipe.roughness);
		double p_0 = task.p_L + (pipe.z_L - pipe.z_0) * task.rho * g +
			(hydraulic_resistance * pipe.length * pow(speed, 2) * task.rho) / (2 * task.pipe.internal_diameter);
		return p_0;
	}
	/// @brief Метод нахождения выходного давления
	/// @return Pвых
	double PQ_task()
	{
		double speed = calc_speed(task.Q, pipe.internal_diameter);
		double Re = calc_Re(speed, pipe.internal_diameter, task.nu);
		double hydraulic_resistance = hydraulic_resistance_isaev(Re, pipe.roughness);
		double p_L = task.p_0 - (pipe.z_L - pipe.z_0) * task.rho * g -
			(hydraulic_resistance * pipe.length * pow(speed, 2) * task.rho) / (2 * task.pipe.internal_diameter);
		return p_L;
	}
	/// @brief Метод нахождения расхода
	/// @return Расход Q
	double PP_task()
	{
		double hydraulic_resistance_prev, speed, Re;
		double hydraulic_resistance = 0.02;
		double hydraulic_resistance_v2 = (pipe.internal_diameter * 2 * g * (((task.p_0 - task.p_L) / (task.rho * g)) + pipe.z_0 - pipe.z_L)) / pipe.length;
		do {
			hydraulic_resistance_prev = hydraulic_resistance;
			speed = sqrt(hydraulic_resistance_v2 / hydraulic_resistance_prev);
			Re = calc_Re(speed, pipe.internal_diameter, task.nu);
			hydraulic_resistance = hydraulic_resistance_isaev(Re, pipe.roughness);
		} while (abs(hydraulic_resistance - hydraulic_resistance_prev) > 0.0002);
		double Q = pi * pow(pipe.internal_diameter, 2) * speed / 4;
		return Q;
	}
};
/// @brief Класс, решающий задачи PQ, QP методом Эйлера
class euler_solver
{
	const my_pipe_parameters& pipe;
	const my_task_parameters& task;
public:
	euler_solver(const my_pipe_parameters& pipe, const my_task_parameters& task) :
		pipe(pipe), task(task)
	{
	}
	/// @brief Метод нахождения входного давления
	/// @return Pвх
	vector<double> euler_solver_PQ()
	{
		double delta_z = (pipe.z_L - pipe.z_0) / (pipe.n - 1);
		double speed = calc_speed(task.Q, pipe.internal_diameter);
		double Re = calc_Re(speed, pipe.internal_diameter, task.nu);
		double hydraulic_resistance = hydraulic_resistance_isaev(Re, pipe.roughness);
		double tau = calc_tau(hydraulic_resistance, task.rho, speed);
		vector<double> p_profile = vector<double>(pipe.n);
		p_profile[0] = task.p_0;
		for (int i = 1; i < pipe.n; i++)
			p_profile[i] = p_profile[i - 1] + pipe.h * ((-4 / pipe.internal_diameter) * tau - task.rho * g * (delta_z / pipe.h));
		return p_profile;
	}
	/// @brief Метод нахождения выходного давления
	/// @return Pвых
	vector<double> euler_solver_QP()
	{
		double delta_z = (pipe.z_L - pipe.z_0) / (pipe.n - 1);
		double speed = calc_speed(task.Q, pipe.internal_diameter);
		double Re = calc_Re(speed, pipe.internal_diameter, task.nu);
		double hydraulic_resistance = hydraulic_resistance_isaev(Re, pipe.roughness);
		double tau = calc_tau(hydraulic_resistance, task.rho, speed);
		vector<double> p_profile = vector<double>(pipe.n);
		p_profile[pipe.n - 1] = task.p_L;
		for (int i = pipe.n - 2; i >= 0; i--)
			p_profile[i] = p_profile[i + 1] - pipe.h * ((-4 / pipe.internal_diameter) * tau - task.rho * g * (delta_z / pipe.h));
		return p_profile;
	}
};
/// @brief Класс, решающий задачу PP методом Ньютона
class newton_solver_PP : public fixed_system_t<1>
{
	const my_pipe_parameters& pipe;
	const my_task_parameters& task;
	using fixed_system_t<1>::var_type;
public:
	newton_solver_PP(const my_pipe_parameters& pipe, const my_task_parameters& task) :
		pipe(pipe), task(task)
	{
	}
	/// @brief Задание функции невязок
	/// @param x Искомая скорость
	/// @return Функция невязок
	var_type residuals(const var_type& x)
	{
		double speed = x;
		double Re = calc_Re(speed, pipe.internal_diameter, task.nu);
		double hydraulic_resistance = hydraulic_resistance_isaev(Re, pipe.roughness);
		return (hydraulic_resistance * (pipe.length * pow(speed, 2) / (pipe.internal_diameter * 2 * g)) + task.p_L / (task.rho * g) + pipe.z_L - task.p_0 / (task.rho * g) - pipe.z_0);
	}
};
/// @brief Класс, решающий задачу PP методом Ньютона поверх Эйлера
class newton_solver_PP_with_euler : public fixed_system_t<1>
{
	const my_pipe_parameters& pipe;
	const my_task_parameters& task;
	using fixed_system_t<1>::var_type;
public:
	newton_solver_PP_with_euler(const my_pipe_parameters& pipe, const my_task_parameters& task) :
		pipe(pipe), task(task)
	{
	}
	/// @brief Задание функции невязок
	/// @param x Искомый расход
	/// @return Функция невязок
	var_type residuals(const var_type& x)
	{
		my_task_parameters temp_task = task; // Временная структура
		temp_task.Q = x; // во временной структуре используем Q для нашего уравнения невязки, эта Q будет идти в солвер
		euler_solver e_solver(pipe, temp_task); // Объявляем переменную класса солвера Эйлером
		vector<double> p_profile = e_solver.euler_solver_QP(); // Считаем профиль давлений Эйлером
		return (p_profile[0] - task.p_0);
	}
};
/// @brief Решение задачи PQ методом простых итераций
TEST(MOC_Solver, Task_1)
{
	simple_pipe_properties simple_pipe;
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 80e3;
	simple_pipe.dx = 1e3;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 100, z_L = 50,
		rho = 870, nu = 15e-6, p_0 = 0.0, p_L = 0.6e6, Q = 0.972;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q };
	simple_iterations_solver simple(pipe,task);
	p_0 = simple.QP_task();
	cout << p_0 << endl;
}
/// @brief Решение задачи PP методом простых итераций
TEST(MOC_Solver, Task_2)
{
	simple_pipe_properties simple_pipe;
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 80e3;
	simple_pipe.dx = 1e3;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 5e6, p_L = 0.8e6, Q = 0.0;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q };
	simple_iterations_solver simple(pipe, task);
	Q = simple.PP_task() * 3600;
	cout << Q << endl;
}
/// @brief Решение задачи QP методом Эйлера
TEST(MOC_Solver, Task_3)
{
	simple_pipe_properties simple_pipe;
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 80e3;
	simple_pipe.dx = 1e3;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 0.0, p_L = 0.6e6, Q = 0.972;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q };
	euler_solver e_solver(pipe, task);
	vector<double> p_profile = e_solver.euler_solver_QP();
	p_0 = p_profile[0];
	cout << p_0 << endl;
}
/// @brief Решение задачи PP методом Ньютона
TEST(MOC_Solver, Task_4)
{
	simple_pipe_properties simple_pipe; 
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 80e3;
	simple_pipe.dx = 1e3;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 5e6, p_L = 0.8e6, Q = 0.0;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q };
	newton_solver_PP n_solver(pipe, task);
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	// Задаем начальное приближение по скорости (м/с)
	double v_approx = 20;
	fixed_newton_raphson<1>::solve_dense(n_solver, { v_approx }, parameters, &result);
	cout << result.argument * pi * pow(pipe.internal_diameter, 2) / 4 * 3600 << endl;
}
/// @brief Решение задачи PP методом Ньютона поверх Эйлера
TEST(MOC_Solver, Task_5)
{
	simple_pipe_properties simple_pipe;
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 80e3;
	simple_pipe.dx = 1e3;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 5e6, p_L = 0.8e6, Q = 0.0;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	my_task_parameters task{ pipe, rho , nu, p_0, p_L, Q};
	newton_solver_PP_with_euler n_solver(pipe, task);
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	// Задаем начальное приближение по расходу (м3/с)
	double Q_approx = 0.5;
	fixed_newton_raphson<1>::solve_dense(n_solver, { Q_approx }, parameters, &result);
	cout << result.argument * 3600 << endl;
}