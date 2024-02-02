clc
clear

% Загрузка данных из файла CSV
filename = 'research_out/diff_p_profile.csv';
data = readtable(filename);

data = dlmread(filename, ';', 0, 0);
km = 0:0.1:100;
minValue = 1.1*min(data(:, 2:end-1), [], 'all');
maxValue = 1.1*max(data(:, 2:end-1), [], 'all');

filename2 = 'research_out/p_profile.csv';
data2 = readtable(filename2);
data2 = dlmread(filename2, ';', 0, 0);
minValue2 = min(data2(:, 2:end-1), [], 'all')-0.1e6;
maxValue2 = max(data2(:, 2:end-1), [], 'all')+0.1e6;


filename3 = 'research_out/nu_profile.csv';
data3 = readtable(filename3);
data3 = dlmread(filename3, ';', 0, 0);
minValue3 = min(data3(:, 2:end-1), [], 'all') - 0.1e-5;
maxValue3 = max(data3(:, 2:end-1), [], 'all') + 0.1e-5;;

% Отображение данных перед началом цикла
figure;
subplot(3, 1, 1);
plot(km, data(1, 2:end-1), 'Color', 'b');
hold on;
xlabel('Труба, км');
ylabel('Разница давлений, Па');
title(['t = ' num2str(data(1,1)) ', с']);
newXLimit = [0, 100];
xlim(newXLimit);
newYLimit = [minValue, maxValue];
ylim(newYLimit);

subplot(3, 1, 2);
plot(km, data2(1, 2:end-1), 'Color', 'b');
xlabel('Труба, км');
ylabel('Давление, Па');
title(['t = ' num2str(data2(1,1)) ', с']);
newXLimit = [0, 100];
xlim(newXLimit);
newYLimit = [minValue2, maxValue2];
ylim(newYLimit);

subplot(3, 1, 3);
plot(km, data3(1, 2:end-1), 'Color', 'b');
xlabel('Труба, км');
ylabel('Вязкость, сСт');
title(['t = ' num2str(data3(1,1)) ', с']);
newXLimit = [0, 100];
xlim(newXLimit);
newYLimit = [minValue3, maxValue3];
ylim(newYLimit);

% Получение кадра для первого кадра
frame = getframe(gcf);
[im, map] = rgb2ind(frame.cdata, 256, 'nodither');
im(1, 1, 1, size(data, 1)) = 0;

% Анимация
for i = 2:size(data(:,1))
    % Отображение данных
    subplot(3, 1, 1);
    plot(km, data(i, 2:end-1), 'Color', 'r');
    hold on;
    plot(km, data(1, 2:end-1), 'Color', 'b');
    newXLimit = [0, 100];
    xlim(newXLimit);
    newYLimit = [minValue, maxValue];
    ylim(newYLimit);
    xlabel('Труба, км');
    ylabel('Разница давлений, Па');
    title(['t = ' num2str(data(i,1)) ', с']);
    hold off
    
    subplot(3, 1, 2);
    plot(km, data2(i, 2:end-1), 'Color', 'b');
    xlabel('Труба, км');
    ylabel('Давление, Па');
    title(['t = ' num2str(data2(i,1)) ', с']);
    newXLimit = [0, 100];
    xlim(newXLimit);
    newYLimit = [minValue2, maxValue2];
    ylim(newYLimit);

    subplot(3, 1, 3);
    plot(km, data3(i, 2:end-1), 'Color', 'b');
    xlabel('Труба, км');
    ylabel('Вязкость, сСт');
    title(['t = ' num2str(data3(i,1)) ', с']);
    newXLimit = [0, 100];
    xlim(newXLimit);
    newYLimit = [minValue3, maxValue3];
    ylim(newYLimit);

    % Получение кадра
    frame = getframe(gcf);
    % Добавляем кадр к гифке
    im(:, :, 1, i) = rgb2ind(frame.cdata, map, 'nodither');

end

% Сохранение гифки в файл
filename = 'отклонение+вязкость.gif';
imwrite(im, map, filename, 'DelayTime', 0.02, 'LoopCount', inf);
disp(['Гифка сохранена в файл: ' filename]);