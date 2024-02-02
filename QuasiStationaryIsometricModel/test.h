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
	
	vector<double> rho_profile = vector<double>(pipe.n, rho);
	vector<double> nu_profile = vector<double>(pipe.n, nu);
	my_task_parameters(my_pipe_parameters& pipe, double rho, double nu, double p_0, double p_L, double Q, vector<double> rho_profile, vector<double> nu_profile) :
		pipe{ pipe }, rho{ rho }, nu{ nu }, p_0{ p_0 }, p_L{ p_L }, Q{ Q }, rho_profile{ rho_profile }, nu_profile{ nu_profile }
	{	
	}
};

struct print_data {
	vector<double> time;           // время с
	vector<double> density;        // вытесняющая плотность кг/м3
	vector<double> viscosity;      // вытесняющая вязкость Ст
	vector<double> inputPressure;  // давление на входе Па
	vector<double> outputPressure; // давление на выходе Па
	vector<double> flowRate;       // расход м3/с
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

/// @brief Функция расчета расхода из скорости
/// @param speed Скорость, (м/с)
/// @param internal_diameter Внутренний диаметр, (м) 
/// @return Расход, (м^3/с)
double calc_flow(double speed, double internal_diameter)
{
	double flow = (pow(internal_diameter, 2) * pi * speed) / 4;
	return flow;
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
/// @brief Класс, решающий задачи PQ, QP методом Эйлера
class euler_solver_with_MOC
{
	const my_pipe_parameters& pipe;
	const my_task_parameters& task;
public:
	euler_solver_with_MOC(const my_pipe_parameters& pipe, const my_task_parameters& task) :
		pipe(pipe), task(task)
	{
	}
	/// @brief Метод нахождения входного давления
	/// @return Pвх
	vector<double> euler_solver_PQ()
	{
		double delta_z = (pipe.z_L - pipe.z_0) / (pipe.n - 1);
		double speed = calc_speed(task.Q, pipe.internal_diameter);
		vector<double> p_profile = vector<double>(pipe.n);
		p_profile[0] = task.p_0;
		for (int i = 1; i < pipe.n; i++)
		{
			double Re = calc_Re(speed, pipe.internal_diameter, task.nu_profile[i]);
			double hydraulic_resistance = hydraulic_resistance_isaev(Re, pipe.roughness);
			double tau = calc_tau(hydraulic_resistance, task.rho_profile[i], speed);
			p_profile[i] = p_profile[i - 1] + pipe.h * ((-4 / pipe.internal_diameter) * tau - task.rho_profile[i] * g * (delta_z / pipe.h));

		}
		return p_profile;
	}
	/// @brief Метод нахождения выходного давления
	/// @return Pвых
	vector<double> euler_solver_QP()
	{
		double delta_z = (pipe.z_L - pipe.z_0) / (pipe.n - 1);
		double speed = calc_speed(task.Q, pipe.internal_diameter);
		vector<double> p_profile = vector<double>(pipe.n);
		p_profile[pipe.n - 1] = task.p_L;
		for (int i = pipe.n - 2; i >= 0; i--)
		{
			double Re = calc_Re(speed, pipe.internal_diameter, task.nu_profile[i]);
			double hydraulic_resistance = hydraulic_resistance_isaev(Re, pipe.roughness);
			double tau = calc_tau(hydraulic_resistance, task.rho_profile[i], speed);
			p_profile[i] = p_profile[i + 1] - pipe.h * ((-4 / pipe.internal_diameter) * tau - task.rho_profile[i] * g * (delta_z / pipe.h));

		}
		return p_profile;
	}
};

class newton_solver_PP_with_euler_with_MOC : public fixed_system_t<1>
{
	const my_pipe_parameters& pipe;
	const my_task_parameters& task;
	using fixed_system_t<1>::var_type;
public:
	newton_solver_PP_with_euler_with_MOC(const my_pipe_parameters& pipe, const my_task_parameters& task) :
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
		euler_solver_with_MOC e_solver(pipe, temp_task); // Объявляем переменную класса солвера Эйлером

		p_profile = e_solver.euler_solver_PQ(); // Считаем профиль давлений Эйлером
		return (p_profile.back() - task.p_L);
	}
	vector<double> get_p_profile() const {
		return p_profile;
	}
private:
	vector<double> p_profile;
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
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, {}, {}};
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
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, {}, {}};
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
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, {}, {}};
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
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, {}, {}};
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
	my_task_parameters task{ pipe, rho , nu, p_0, p_L, Q, {}, {}};
	newton_solver_PP_with_euler n_solver(pipe, task);
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	// Задаем начальное приближение по расходу (м3/с)
	double Q_approx = 0.5;
	fixed_newton_raphson<1>::solve_dense(n_solver, { Q_approx }, parameters, &result);
	cout << result.argument * 3600 << endl;
}

/// @brief тест не рабочий buffer.advance
TEST(MOC_Solver, Task_6)
{
	simple_pipe_properties simple_pipe;
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 500;
	simple_pipe.dx = 10;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 100, z_L = 50,
		rho = 900, nu = 15e-6, p_0 = 6e6, p_L = 0.0, speed = 0.5;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	double Q = calc_flow(speed, pipe.internal_diameter);
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, {}, {}};

	double rho_in = 800;
	double rho_out = 900;
	double nu_in = 10e-6;
	double nu_out = 15e-6;
	double T = 200;

	pipe_properties_t simple_pipe_properties = pipe_properties_t::build_simple_pipe(simple_pipe);

	string path = "./";
	
	typedef composite_layer_t<profile_collection_t<3>,
		moc_solver<1>::specific_layer> single_var_moc_t;

	vector<double> Q_profile(simple_pipe_properties.profile.getPointCount(), Q);
	PipeQAdvection advection_model(simple_pipe_properties, Q_profile);

	const auto& x = advection_model.get_grid();
	double dx = x[1] - x[0];
	double v = advection_model.getEquationsCoeffs(0, 0);
	double dt_ideal = abs(dx / v);

	double Cr = 0.5;

	ring_buffer_t<single_var_moc_t> buffer(2, simple_pipe_properties.profile.getPointCount());
	buffer.advance(+1);
	single_var_moc_t& prev = buffer.previous();
	single_var_moc_t& next = buffer.current();

	euler_solver_with_MOC e_solver(pipe, task);

	auto& rho_initial = prev.vars.point_double[0];
	rho_initial = vector<double>(rho_initial.size(), rho);
	task.rho_profile = rho_initial;

	auto& nu_initial = prev.vars.point_double[1];
	nu_initial = vector<double>(nu_initial.size(), nu);
	task.nu_profile = nu_initial;

	prev.vars.point_double[2] = e_solver.euler_solver_PQ();
	
	double t = 0; // текущее время
	double dt = Cr * dt_ideal; // время в долях от Куранта
	
	std::stringstream filename;
	filename << path << "output.csv";
	std::ofstream output(filename.str());

	size_t N = static_cast<int>(T / dt);
	for (size_t index = 0; index < N; ++index) {

		if (index == 0) {
			single_var_moc_t& prev = buffer.previous();
			prev.vars.print(t, output);
		}
		t += dt;

		moc_layer_wrapper<1> moc_prev_rho(prev.vars.point_double[0],
			std::get<0>(prev.specific));
		moc_layer_wrapper<1> moc_next_rho(next.vars.point_double[0],
			std::get<0>(next.specific));

		moc_solver<1> solver(advection_model, moc_prev_rho, moc_next_rho);

		solver.step_optional_boundaries(dt, rho_in, rho_out);	

		moc_layer_wrapper<1> moc_prev_nu(prev.vars.point_double[1],
			std::get<0>(prev.specific));
		moc_layer_wrapper<1> moc_next_nu(next.vars.point_double[1],
			std::get<0>(next.specific));

		moc_solver<1> solver2(advection_model, moc_prev_nu, moc_next_nu);

		solver2.step_optional_boundaries(dt, nu_in, nu_out);
		next = buffer.current();
		task.rho_profile = next.vars.point_double[0];
		task.nu_profile = next.vars.point_double[1];

		next.vars.point_double[2] = e_solver.euler_solver_PQ();
		
		next.vars.print(t, output);


		buffer.advance(+1);

	}
	output.flush();
	output.close();
}


// @brief класс, созданный для решения транспортного уравнения
class simple_moc_solver {
public:
	const my_pipe_parameters& pipe;
	const my_task_parameters& task;
	simple_moc_solver(const my_pipe_parameters& pipe, const my_task_parameters& task, vector<vector<double>>& layer_prev, vector<vector<double>>& layer_curr) :
		pipe(pipe), task(task), layer_prev(layer_prev), layer_curr(layer_curr)
	{
	}

	double prepare_step(double time_step = std::numeric_limits<double>::quiet_NaN()) {
		auto& values = layer_prev;

		double eigen_value = calc_speed(task.Q, pipe.internal_diameter);
	
		double courant_step = pipe.h / eigen_value;
		if (std::isnan(time_step) || time_step > courant_step) {
			time_step = courant_step;
		}
		return time_step;
	}

	/// @brief метод, делающий реальный расчет, в него передаются два граничных условия
	/// @param left_value левое граничное условие
	/// @param right_value правое граничное условие
	void step(double time_step, double dt, vector<vector<double>> left_value, vector<vector<double>> right_value)
	{
// не знаю как и куда всунуть расчитанный курант, но тут он всегда 1, поэтому норм
		double speed = calc_speed(task.Q, pipe.internal_diameter);
		double dt_ideal = abs(pipe.h / speed);
		double Cr = speed * dt / pipe.h;
// основа (для куранта 1)
		for (size_t num_prof = 0; num_prof < layer_prev.size(); num_prof++)
		{
			if (task.Q > 0)
			{
				for (int l = 0; l < pipe.n - 1; l++)
					layer_curr[num_prof][l + 1] = layer_prev[num_prof][l];
				layer_curr[num_prof][0] = left_value[num_prof][time_step-1];
			}
			else
			{
				for (int l = pipe.n - 1; l > 0; l--)
					layer_curr[num_prof][l - 1] = layer_prev[num_prof][l];
				layer_curr[num_prof][pipe.n] = right_value[num_prof][time_step-1];
			}
		}

	}
private:
	vector<vector<double>>& layer_prev;
	vector<vector<double>>& layer_curr;
};

double linear_interpolator(vector<double> original_time, vector<double> original_value, double new_time_step) 
{
	size_t index1 = 0;
	size_t index2 = 1;
	while (original_time[index2] < new_time_step) 
	{
		++index1;
		++index2;
	}
	double t1 = original_time[index1];
	double t2 = original_time[index2];
	double value1 = original_value[index1];
	double value2 = original_value[index2];
	return value1 + (value2 - value1) * (new_time_step - t1) / (t2 - t1);
}

void print_data_to_csv(const vector<double>& time,
	const vector<double>& rho_and_nu_in_0,
	const vector<double>& rho_and_nu_in_1,
	const vector<double>& time_p_in,
	const vector<double>& time_p_out,
	const vector<double>& time_Q,
	const wstring& filename)
{

	// Определяем максимальную длину вектора
	size_t maxLength = max({ time.size(), rho_and_nu_in_0.size(), rho_and_nu_in_1.size(),
									  time_p_in.size(), time_p_out.size(), time_Q.size() });

	// Открываем файл для записи
	ofstream file(filename);

	// Проверяем, открыт ли файл успешно
	if (file.is_open()) {
		// Записываем заголовки столбцов
		file << "Time; Density; Viscosity; Time Pressure In; Time Pressure Out; Time Flow Rate\n";

		// Записываем данные из векторов
		for (size_t i = 0; i < maxLength; ++i) {
			// Если индекс находится в пределах длины вектора, записываем значение, иначе записываем пустую ячейку
			file << (i < time.size() ? to_string(time[i]) : "") << ";"
				<< (i < rho_and_nu_in_0.size() ? to_string(rho_and_nu_in_0[i]) : "") << ";"
				<< (i < rho_and_nu_in_1.size() ? to_string(rho_and_nu_in_1[i]) : "") << ";"
				<< (i < time_p_in.size() ? to_string(time_p_in[i]) : "") << ";"
				<< (i < time_p_out.size() ? to_string(time_p_out[i]) : "") << ";"
				<< (i < time_Q.size() ? to_string(time_Q[i]) : "") << "\n";
		}
		// Закрываем файл
		file.close();
	}
}


void print_layers(const double dt, 
	const vector<double>& layer, 
	const wstring& filename)
{
	ofstream  file(filename, ios::app);
	if (dt == 0)
	{
		// Если файл существует, очищаем его содержимое
		if (file.is_open()) {
			file.close();
			file.open(filename, ios::out | ios::trunc);
		}
		// Если файл не существует, создаем новый
		else {
			file.open(filename, ios::out);
		}
		// Файл существует, но file.is_open() возвращает False
		if (!file.is_open()) {
			file.clear();  // Очищаем флаг ошибки
			file.open(filename, ios::out | ios::trunc);
		}
	}
	if (file.is_open()) {
		file << to_string(dt) << ";";
		for (int j = 0; j < layer.size(); j++)
		{
			file << to_string(layer[j]) << ";";
		}
		file << "\n";
		file.close();
	}
}

TEST(Block_3, Task_2)
{
	simple_pipe_properties simple_pipe;
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 500;
	simple_pipe.dx = 10;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 900, nu = 15e-6, p_0 = 6e6, p_L = 0.0, speed = 0.5;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	double Q = calc_flow(speed, pipe.internal_diameter);
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, {}, {}};
	
	std::srand(std::time(nullptr));
	vector<double> time_row;
	for (double time = 0.0; time <= 360.0; time += 10.0)
		time_row.push_back(time);
	vector<double> time_rho_in_row = vector<double>(time_row.size(), rho);
	for (size_t i = 5; i <= 36; i++)
	{
		time_rho_in_row[i] = rho - abs(rho * 0.05 * std::rand() / RAND_MAX);
	}
	vector<double> time_nu_in_row = vector<double>(time_row.size(), nu);	
	for (size_t i = 5; i <= 36; i++)
	{
		time_nu_in_row[i] = nu - abs(nu * 0.2 * std::rand() / RAND_MAX);
	}
	vector<double> time_rho_out_row = vector<double>(time_row.size(), rho);
	vector<double> time_nu_out_row = vector<double>(time_row.size(), nu);
	vector<double> time_p_in_row = vector<double>(time_row.size(), p_0);
	for (size_t i = 5; i <= 36; i++)
	{
		time_p_in_row[i] = p_0 - abs(p_0 * 2e-2 * std::rand() / RAND_MAX);
	}
	vector<double> time_Q_row = vector<double>(time_row.size(), Q);
	for (size_t i = 5; i <= 36; i++)
	{
		time_Q_row[i] = Q - abs(Q * 0.1 * std::rand() / RAND_MAX);
	}

	vector<vector<double>> layer = vector<vector<double>>(2, vector<double>(pipe.n));
	ring_buffer_t<vector<vector<double>>> buffer(2, layer);

	buffer.advance(+1);
	buffer.previous()[0] = vector<double>(buffer.previous()[0].size(), rho);
	buffer.previous()[1] = vector<double>(buffer.previous()[1].size(), nu);

	vector<double> new_time_row, new_time_p_in_row, new_time_p_out_row, new_time_Q_row, p_profile, initial_p_profile;
	vector<double> diff_p_profile = vector<double>(pipe.n);
	vector<vector<double>> rho_and_nu_in = vector<vector<double>>(2);
	vector<vector<double>> rho_and_nu_out = vector<vector<double>>(2);

	double dt = 0;

	euler_solver_with_MOC e_solver(pipe, task);

	wstring folder_path = L"research\\2024_02_block_3\\task_2\\research_out";
	wstring p_profile_file = folder_path + L"\\p_profile.csv";
	wstring rho_profile_file = folder_path + L"\\rho_profile.csv";
	wstring nu_profile_file = folder_path + L"\\nu_profile.csv";
	wstring diff_p_profile_file = folder_path + L"\\diff_p_profile.csv";

	do {
		new_time_row.push_back(dt);
		task.Q = linear_interpolator(time_row, time_Q_row, dt);
		new_time_Q_row.push_back(task.Q);
		rho_and_nu_in[0].push_back(linear_interpolator(time_row, time_rho_in_row, dt));
		rho_and_nu_in[1].push_back(linear_interpolator(time_row, time_nu_in_row, dt));

		rho_and_nu_out[0].push_back(linear_interpolator(time_row, time_rho_out_row, dt));
		rho_and_nu_out[1].push_back(linear_interpolator(time_row, time_nu_out_row, dt));

		new_time_p_in_row.push_back(linear_interpolator(time_row, time_p_in_row, dt));

		simple_moc_solver simple_moc(pipe, task, buffer.previous(), buffer.current());

		simple_moc.step(new_time_row.size(), simple_moc.prepare_step(), rho_and_nu_in, rho_and_nu_out);
		task.rho_profile = buffer.current()[0];
		task.nu_profile = buffer.current()[1];
		task.p_0 = new_time_p_in_row.back();
		task.Q = new_time_Q_row.back();
		p_profile = e_solver.euler_solver_PQ();
		if (dt == 0)
			initial_p_profile = p_profile;
		new_time_p_out_row.push_back(p_profile.back());
		std::transform(initial_p_profile.begin(), initial_p_profile.end(), p_profile.begin(), diff_p_profile.begin(),
			[](double initial, double current) {return initial - current;  });
		print_layers(dt, p_profile, p_profile_file);
		print_layers(dt, buffer.current()[0], rho_profile_file);
		print_layers(dt, buffer.current()[1], nu_profile_file);
		print_layers(dt, diff_p_profile, diff_p_profile_file);
		buffer.advance(+1);
		dt += simple_moc.prepare_step(); // здесь используется task.Q для расчета шага
		//все готово для следующего шага, интерполированная скорость есть
	} while (dt < time_row.back());	
	wstring filename_initial = folder_path + L"\\initial_data.csv";
	print_data_to_csv(
		time_row,
		time_rho_in_row,
		time_nu_in_row,
		time_p_in_row,
		{},
		time_Q_row,
		filename_initial
	);	
	wstring filename_final = folder_path + L"\\final_data.csv";
	print_data_to_csv(
		new_time_row,
		rho_and_nu_in[0],
		rho_and_nu_in[1],
		new_time_p_in_row,
		new_time_p_out_row,
		new_time_Q_row,
		filename_final
	);
}


TEST(Block_3, Task_3)
{
	simple_pipe_properties simple_pipe;
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 500;
	simple_pipe.dx = 10;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 900, nu = 15e-6, p_0 = 6e6, p_L = 5.557e6, Q = 0;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, {}, {} };

	vector<double> time_row = { 0, 60, 120, 180, 240, 300, 360 };
	vector<double> time_rho_in_row = { rho, 880, 880, 890, 890, 880, 880 };
	vector<double> time_nu_in_row = { nu, 13e-6, 13e-6, 14e-6, 14e-6, 13e-6, 13e-6 };

	vector<double> time_rho_out_row = { rho, 880, 880, 890, 890, 880, 880 };
	vector<double> time_nu_out_row = { nu, 13e-6, 13e-6, 14e-6, 14e-6, 13e-6, 13e-6 };

	vector<double> time_p_in_row = { p_0, 5.8e6, 5.8e6, 5.9e6, 5.9e6, 5.8e6, 5.8e6 };
	vector<double> time_p_out_row = { p_L, 5.357e6, 5.357e6, 5.458e6, 5.458e6, 5.359e6, 5.359e6 };

	vector<vector<double>> layer = vector<vector<double>>(2, vector<double>(pipe.n));
	ring_buffer_t<vector<vector<double>>> buffer(2, layer);

	buffer.advance(+1);
	buffer.previous()[0] = vector<double>(buffer.previous()[0].size(), rho);
	buffer.previous()[1] = vector<double>(buffer.previous()[1].size(), nu);
	
	task.rho_profile = buffer.previous()[0];
	task.nu_profile = buffer.previous()[1];

	vector<double> new_time_row, new_time_p_in_row, new_time_p_out_row, new_time_Q_row, p_profile;
	vector<vector<double>> rho_and_nu_in = vector<vector<double>>(2);
	vector<vector<double>> rho_and_nu_out = vector<vector<double>>(2);

	double dt = 0;
	
	newton_solver_PP_with_euler_with_MOC n_solver(pipe, task);
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;

	double Q_approx = 0.19;
	wstring folder_path = L"research\\2024_02_block_3\\task_3\\research_out";
	wstring p_profile_file = folder_path + L"\\p_profile.csv";
	wstring rho_profile_file = folder_path + L"\\rho_profile.csv";
	wstring nu_profile_file = folder_path + L"\\nu_profile.csv";
	do {
		new_time_row.push_back(dt);
		rho_and_nu_in[0].push_back(linear_interpolator(time_row, time_rho_in_row, dt));
		rho_and_nu_in[1].push_back(linear_interpolator(time_row, time_nu_in_row, dt));
		rho_and_nu_out[0].push_back(linear_interpolator(time_row, time_rho_out_row, dt));
		rho_and_nu_out[1].push_back(linear_interpolator(time_row, time_nu_out_row, dt));

		new_time_p_in_row.push_back(linear_interpolator(time_row, time_p_in_row, dt));
		new_time_p_out_row.push_back(linear_interpolator(time_row, time_p_out_row, dt));

		fixed_newton_raphson<1>::solve_dense(n_solver, { Q_approx }, parameters, &result);
		task.Q = result.argument;
		new_time_Q_row.push_back(task.Q);
		simple_moc_solver simple_moc(pipe, task, buffer.previous(), buffer.current());		
		simple_moc.step(new_time_row.size(), simple_moc.prepare_step(), rho_and_nu_in, rho_and_nu_out);
		task.rho_profile = buffer.current()[0];
		task.nu_profile = buffer.current()[1];
		task.rho = task.rho_profile[0];
		task.nu = task.nu_profile[0];
		task.p_0 = new_time_p_in_row.back();
		task.p_L = new_time_p_out_row.back();		
		p_profile = n_solver.get_p_profile();		
		print_layers(dt, p_profile, p_profile_file);
		print_layers(dt, buffer.current()[0], rho_profile_file);
		print_layers(dt, buffer.current()[1], nu_profile_file);
		buffer.advance(+1);
		dt += simple_moc.prepare_step(); // здесь используется task.Q для расчета шага
		//все готово для следующего шага, интерполированная скорость есть
	} while (dt < time_row.back());
	wstring filename_initial = folder_path + L"\\initial_data.csv";
	print_data_to_csv(
		time_row,
		time_rho_in_row,
		time_nu_in_row,
		time_p_in_row,
		time_p_out_row,
		{},
		filename_initial
	);
	wstring filename_final = folder_path + L"\\final_data.csv";
	print_data_to_csv(
		new_time_row,
		rho_and_nu_in[0],
		rho_and_nu_in[1],
		new_time_p_in_row,
		new_time_p_out_row,
		new_time_Q_row,
		filename_final
	);
}

TEST(Block_3, Task_2_One_Parameter)
{
	simple_pipe_properties simple_pipe;
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 100e3;
	simple_pipe.dx = 100;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 900, nu = 15e-6, p_0 = 6e6, p_L = 0.0, speed = 0.5;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	double Q = calc_flow(speed, pipe.internal_diameter);
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, {}, {} };
	
	vector<double> time_row;
	for (double time = 0.0; time <= 200000.0; time += 100.0) 
		time_row.push_back(time);
	
	vector<double> time_rho_in_row = vector<double>(time_row.size(), rho);
	vector<double> time_nu_in_row = vector<double>(time_row.size(), nu);
	std::srand(std::time(nullptr));
	for (size_t i = 10; i <= 2000; i++)
	{
		time_nu_in_row[i] = nu  + nu * 0.2 * std::rand() / RAND_MAX;
	}
	vector<double> time_rho_out_row = vector<double>(time_row.size(), rho);
	vector<double> time_nu_out_row = vector<double>(time_row.size(), nu);

	vector<double> time_p_in_row = vector<double>(time_row.size(), p_0);
	vector<double> time_Q_row = vector<double>(time_row.size(), Q);


	vector<vector<double>> layer = vector<vector<double>>(2, vector<double>(pipe.n));
	ring_buffer_t<vector<vector<double>>> buffer(2, layer);

	buffer.advance(+1);
	buffer.previous()[0] = vector<double>(buffer.previous()[0].size(), rho);
	buffer.previous()[1] = vector<double>(buffer.previous()[1].size(), nu);

	vector<double> new_time_row, new_time_p_in_row, new_time_p_out_row, new_time_Q_row, p_profile, initial_p_profile;
	vector<double> diff_p_profile = vector<double>(pipe.n);
	vector<vector<double>> rho_and_nu_in = vector<vector<double>>(2);
	vector<vector<double>> rho_and_nu_out = vector<vector<double>>(2);

	double dt = 0;

	euler_solver_with_MOC e_solver(pipe, task);

	wstring folder_path = L"research\\2024_02_block_3\\task_2_change_nu\\research_out";
	wstring p_profile_file = folder_path + L"\\p_profile.csv";
	wstring rho_profile_file = folder_path + L"\\rho_profile.csv";
	wstring nu_profile_file = folder_path + L"\\nu_profile.csv";
	wstring diff_p_profile_file = folder_path + L"\\diff_p_profile.csv";

	do {
		new_time_row.push_back(dt);
		task.Q = linear_interpolator(time_row, time_Q_row, dt);
		new_time_Q_row.push_back(task.Q);
		rho_and_nu_in[0].push_back(linear_interpolator(time_row, time_rho_in_row, dt));
		rho_and_nu_in[1].push_back(linear_interpolator(time_row, time_nu_in_row, dt));

		rho_and_nu_out[0].push_back(linear_interpolator(time_row, time_rho_out_row, dt));
		rho_and_nu_out[1].push_back(linear_interpolator(time_row, time_nu_out_row, dt));

		new_time_p_in_row.push_back(linear_interpolator(time_row, time_p_in_row, dt));

		simple_moc_solver simple_moc(pipe, task, buffer.previous(), buffer.current());

		simple_moc.step(new_time_row.size(), simple_moc.prepare_step(), rho_and_nu_in, rho_and_nu_out);
		task.rho_profile = buffer.current()[0];
		task.nu_profile = buffer.current()[1];
		task.p_0 = new_time_p_in_row.back();
		task.Q = new_time_Q_row.back();
		p_profile = e_solver.euler_solver_PQ();
		if (dt == 0)
			initial_p_profile = p_profile;
		new_time_p_out_row.push_back(p_profile.back());
		std::transform(initial_p_profile.begin(), initial_p_profile.end(), p_profile.begin(), diff_p_profile.begin(),
			[](double initial, double current) {return initial - current;  });
		print_layers(dt, p_profile, p_profile_file);
		print_layers(dt, buffer.current()[0], rho_profile_file);
		print_layers(dt, buffer.current()[1], nu_profile_file);
		print_layers(dt, diff_p_profile, diff_p_profile_file);
		buffer.advance(+1);
		dt += simple_moc.prepare_step(); // здесь используется task.Q для расчета шага
		//все готово для следующего шага, интерполированная скорость есть
	} while (dt <= time_row.back());
	wstring filename_initial = folder_path + L"\\initial_data.csv";
	print_data_to_csv(
		time_row,
		time_rho_in_row,
		time_nu_in_row,
		time_p_in_row,
		{},
		time_Q_row,
		filename_initial
	);
	wstring filename_final = folder_path + L"\\final_data.csv";
	print_data_to_csv(
		new_time_row,
		rho_and_nu_in[0],
		rho_and_nu_in[1],
		new_time_p_in_row,
		new_time_p_out_row,
		new_time_Q_row,
		filename_final
	);
}


TEST(Block_3, Task_3_One_Parameter)
{
	simple_pipe_properties simple_pipe;
	simple_pipe.diameter = 0.72;
	simple_pipe.length = 100e3;
	simple_pipe.dx = 100;
	double delta_d = 0.01, abs_roughness = 15e-6, z_0 = 50, z_L = 100,
		rho = 870, nu = 15e-6, p_0 = 6e6, p_L = 5.16e6, Q = 0;
	my_pipe_parameters pipe{ simple_pipe.length, simple_pipe.diameter, delta_d, abs_roughness, z_0, z_L, simple_pipe.dx };
	my_task_parameters task{ pipe, rho, nu, p_0, p_L, Q, {}, {} };
	vector<double> time_row;
	for (double time = 0.0; time <= 200000.0; time += 100.0)
		time_row.push_back(time);
	vector<double> time_rho_in_row = vector<double>(time_row.size(), rho);
	vector<double> time_nu_in_row = vector<double>(time_row.size(), nu);
	std::srand(std::time(nullptr));
	for (size_t i = 10; i <= 500; i++)
	{
		time_rho_in_row[i] = 865;
	}

	for (size_t i = 501; i <= 1000; i++)
	{
		time_rho_in_row[i] = 860;
	}

	for (size_t i = 10; i <= 500; i++)
	{
		time_nu_in_row[i] = 14e-6;
	}

	for (size_t i = 501; i <= 2000; i++)
	{
		time_nu_in_row[i] = 13e-6;
	}
	vector<double> time_rho_out_row = vector<double>(time_row.size(), rho);
	vector<double> time_nu_out_row = vector<double>(time_row.size(), nu);
	vector<double> time_p_in_row = vector<double>(time_row.size(), p_0);
	vector<double> time_p_out_row = vector<double>(time_row.size(), p_L);


	vector<vector<double>> layer = vector<vector<double>>(2, vector<double>(pipe.n));
	ring_buffer_t<vector<vector<double>>> buffer(2, layer);

	buffer.advance(+1);
	buffer.previous()[0] = vector<double>(buffer.previous()[0].size(), rho);
	buffer.previous()[1] = vector<double>(buffer.previous()[1].size(), nu);

	task.rho_profile = buffer.previous()[0];
	task.nu_profile = buffer.previous()[1];

	vector<double> new_time_row, new_time_p_in_row, new_time_p_out_row, new_time_Q_row, p_profile, initial_p_profile;
	vector<double> diff_p_profile = vector<double>(pipe.n);
	vector<vector<double>> rho_and_nu_in = vector<vector<double>>(2);
	vector<vector<double>> rho_and_nu_out = vector<vector<double>>(2);

	double dt = 0;

	newton_solver_PP_with_euler_with_MOC n_solver(pipe, task);
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;

	double Q_approx = 0.19;
	wstring folder_path = L"research\\2024_02_block_3\\task_3_change_rho_and_nu\\research_out";
	wstring p_profile_file = folder_path + L"\\p_profile.csv";
	wstring rho_profile_file = folder_path + L"\\rho_profile.csv";
	wstring nu_profile_file = folder_path + L"\\nu_profile.csv";
	wstring diff_p_profile_file = folder_path + L"\\diff_p_profile.csv";
	do {
		new_time_row.push_back(dt);
		rho_and_nu_in[0].push_back(linear_interpolator(time_row, time_rho_in_row, dt));
		rho_and_nu_in[1].push_back(linear_interpolator(time_row, time_nu_in_row, dt));
		rho_and_nu_out[0].push_back(linear_interpolator(time_row, time_rho_out_row, dt));
		rho_and_nu_out[1].push_back(linear_interpolator(time_row, time_nu_out_row, dt));

		new_time_p_in_row.push_back(linear_interpolator(time_row, time_p_in_row, dt));
		new_time_p_out_row.push_back(linear_interpolator(time_row, time_p_out_row, dt));

		fixed_newton_raphson<1>::solve_dense(n_solver, { Q_approx }, parameters, &result);
		task.Q = result.argument;
		new_time_Q_row.push_back(task.Q);
		simple_moc_solver simple_moc(pipe, task, buffer.previous(), buffer.current());
		simple_moc.step(new_time_row.size(), simple_moc.prepare_step(), rho_and_nu_in, rho_and_nu_out);
		task.rho_profile = buffer.current()[0];
		task.nu_profile = buffer.current()[1];
		task.rho = task.rho_profile[0];
		task.nu = task.nu_profile[0];
		task.p_0 = new_time_p_in_row.back();
		task.p_L = new_time_p_out_row.back();
		p_profile = n_solver.get_p_profile();
		if (dt == 0)
			initial_p_profile = p_profile;
		std::transform(initial_p_profile.begin(), initial_p_profile.end(), p_profile.begin(), diff_p_profile.begin(),
			[](double initial, double current) {return initial - current;  });
		print_layers(dt, p_profile, p_profile_file);
		print_layers(dt, buffer.current()[0], rho_profile_file);
		print_layers(dt, buffer.current()[1], nu_profile_file);
		print_layers(dt, diff_p_profile, diff_p_profile_file);
		buffer.advance(+1);
		dt += simple_moc.prepare_step(); // здесь используется task.Q для расчета шага
		//все готово для следующего шага, интерполированная скорость есть
	} while (dt < time_row.back());
	wstring filename_initial = folder_path + L"\\initial_data.csv";
	print_data_to_csv(
		time_row,
		time_rho_in_row,
		time_nu_in_row,
		time_p_in_row,
		time_p_out_row,
		{},
		filename_initial
	);
	wstring filename_final = folder_path + L"\\final_data.csv";
	print_data_to_csv(
		new_time_row,
		rho_and_nu_in[0],
		rho_and_nu_in[1],
		new_time_p_in_row,
		new_time_p_out_row,
		new_time_Q_row,
		filename_final
	);
}