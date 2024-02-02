% Загрузка данных из файла CSV
filename1 = 'C:\Users\Egor\source\repos\QuasiStationaryIsometricModel\QuasiStationaryIsometricModel\final_data.csv';
filename2 = 'C:\Users\Egor\source\repos\QuasiStationaryIsometricModel\QuasiStationaryIsometricModel\initial_data.csv';
data1 = readtable(filename1);
data2 = readtable(filename2);

% Создание нового графического окна
figure;

% Построение первого графика
subplot(5, 1, 1);
plot(data1.Time, data1.Density, Color='b');
hold on;
plot(data2.Time, data2.Density, Marker="*", Color='r', LineStyle='none');
title('График 1: Временной ряд плотности');
xlabel('Время');
ylabel('Плотность');
legend('Интерполированные данные','Начальные краевые условия')

% Построение второго графика
subplot(5, 1, 2);
plot(data1.Time, data1.Viscosity, Color='b');
hold on;
plot(data2.Time, data2.Viscosity, Marker="*", Color='r', LineStyle='none');
title('График 2: Временной ряд вязкости');
xlabel('Время');
ylabel('Вязкость');
legend('Интерполированные данные','Начальные краевые условия')

% Построение третьего графика
subplot(5, 1, 3);
plot(data1.Time, data1.TimePressureIn, Color='b');
hold on;
plot(data2.Time, data2.TimePressureIn, Marker="*", Color='r', LineStyle='none');
title('График 3: Временной ряд давления на входе');
xlabel('Время');
ylabel('Давление на входе');
legend('Интерполированные данные','Начальные краевые условия')

% Построение четвертого графика
subplot(5, 1, 4);
plot(data1.Time, data1.TimePressureOut, Color='b');
hold on;
plot(data2.Time, data2.TimePressureOut, Marker="*", Color='r', LineStyle='none');
title('График 4: Временной ряд давления на выходе');
xlabel('Время');
ylabel('Давление на выходе');
legend('Интерполированные данные','Начальные краевые условия')

% Построение пятого графика
subplot(5, 1, 5);
plot(data1.Time, data1.TimeFlowRate, Color='b');
hold on;
plot(data2.Time, data2.TimeFlowRate, Marker="*", Color='r', LineStyle='none');
title('График 5: Временной ряд расхода');
xlabel('Время');
ylabel('Расход');
legend('Интерполированные данные','Начальные краевые условия')



% Регулировка размера окна
figure_size = [0, 0, 1920, 1080];

set(gcf, 'Position', figure_size);
%figure_size = [0, 0, 960, 1080];
%set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'output_plot.png');