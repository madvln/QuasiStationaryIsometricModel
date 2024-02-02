% Загрузка данных из файла
filename = 'C:\Users\Egor\source\repos\QuasiStationaryIsometricModel\QuasiStationaryIsometricModel\p_profile.csv';
global data;
data = dlmread(filename, ';', 1, 0);
data = data(:, 2:end-1);
km = 0:10:500;
minValue = min(data, [], 'all')
maxValue = max(data, [], 'all')
% Создание гифки
figure;
for i = 1:size(data, 1)
    % Отображение данных
    
    xlabel('Труба, м');
    ylabel('Давление, Па');
    newXLimit = [0, 500]; % Замените на свой желаемый диапазон
    xlim(newXLimit);
    newYLimit = [minValue, maxValue]; % Замените на свой желаемый диапазон
    ylim(newYLimit);
    % Получение кадра
    frame = getframe(gcf);
    plot(km, data(1, :), '-o',Color='b');
    hold on;
    plot(km, data(i, :), '-o',Color='r');
    
    % Добавление кадра к гифке
    if i == 1
        % Если это первый кадр, создаем гифку
        [im, map] = rgb2ind(frame.cdata, 256, 'nodither');
        im(1, 1, 1, size(data, 1)) = 0;
    else
        % Добавляем кадр к гифке
        im(:, :, 1, i) = rgb2ind(frame.cdata, map, 'nodither');
    end
    hold off;
end

% Сохранение гифки в файл
filename = 'отклонение.gif';
imwrite(im, map, filename, 'DelayTime', 0.5, 'LoopCount', inf);
disp(['Гифка сохранена в файл: ' filename]);