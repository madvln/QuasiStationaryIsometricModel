% Загрузка данных из файла CSV
filename1 = 'research_out/final_data.csv';
filename2 = 'research_out/initial_data.csv';

data1 = readtable(filename1);
data2 = readtable(filename2);

% Создание нового графического окна
figure;

% Построение первого графика
subplot(5, 1, 1);
plot(data2.TimeDensity, data2.Density, Color='0 0 0', LineWidth=2, Marker='.',MarkerSize=20);
hold on;
plot(data1.TimeDensity, data1.Density, Marker=".", Color='1 0.549 0', LineStyle='none', MarkerSize=20);
ylim([890, 910]);
title('График 1: Временной ряд плотности');
xlabel('Время');
ylabel('Плотность');
legend('Исходные краевые условия','Интерполированные данные')

% Построение второго графика
subplot(5, 1, 2);
plot(data2.TimeViscosity, data2.Viscosity, Color='0 0 0', LineWidth=2, Marker='.',MarkerSize=20);
hold on;
plot(data1.TimeViscosity, data1.Viscosity, Marker=".", Color='1 0.549 0', LineStyle='none', MarkerSize=20);
ylim([0.5e-5, 2.5e-5]);
title('График 2: Временной ряд вязкости');
xlabel('Время');
ylabel('Вязкость');
legend('Исходные краевые условия','Интерполированные данные')

% Построение третьего графика
subplot(5, 1, 3);
plot(data2.TimeFlowRate, data2.FlowRate, Color='0 0 0', LineWidth=2, Marker='.',MarkerSize=20);
hold on;
plot(data1.TimeFlowRate, data1.FlowRate, Marker=".", Color='1 0.549 0', LineStyle='none', MarkerSize=20);
ylim([0, 0.4]);
title('График 3: Временной ряд расхода');
xlabel('Время');
ylabel('Расход');
legend('Исходные краевые условия','Интерполированные данные')

% Построение четвертого графика
subplot(5, 1, 4);
plot(data2.TimePressureIn, data2.PressureIn, Color='0 0 0', LineWidth=2, Marker='.',MarkerSize=20);
hold on;
plot(data1.TimePressureIn, data1.PressureIn, Marker=".", Color='1 0.549 0', LineStyle='none', MarkerSize=20);
ylim([5.99e6, 6.01e6]);
title('График 4: Временной ряд давления на входе');
xlabel('Время');
ylabel('Давление на входе');
newXLimit = [0, 400];
legend('Исходные краевые условия','Интерполированные данные')

% Построение пятого графика
subplot(5, 1, 5);
plot(data1.TimePressureOut, data1.PressureOut, Marker=".", Color='1 0.549 0', LineStyle='none', MarkerSize=20);
ylim([5.3e6, 5.7e6]);
title('График 5: Временной ряд давления на выходе');
xlabel('Время');
ylabel('Давление на выходе');
newXLimit = [0, 400];
xlim(newXLimit);
legend('Расчетные значения')


% Регулировка размера окна
figure_size = [0, 0, 1920, 1080];

set(gcf, 'Position', figure_size);
%figure_size = [0, 0, 960, 1080];
%set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'output_plot.png');