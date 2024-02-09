function main()
    clc;
    clear;

    % Загрузка данных из файлов CSV
    [data, km] = loadData('research_out/diff_p_profile.csv');
    [data2, ~] = loadData('research_out/p_profile.csv');
    [data3, ~] = loadData('research_out/rho_profile.csv');
    
    minValue = min(data(:, 2:end-1), [], 'all') - max(data(:, 2:end-1), [], 'all')*0.1;
    maxValue = max(data(:, 2:end-1), [], 'all') + max(data(:, 2:end-1), [], 'all')*0.1;
    minValue2 = min(data2(:, 2:end-1), [], 'all')-0.1e6;
    maxValue2 = max(data2(:, 2:end-1), [], 'all')+0.1e6;
    minValue3 = min(data3(:, 2:end-1), [], 'all')-10;
    maxValue3 = max(data3(:, 2:end-1), [], 'all')+10;
    minValue4 = min(data2(:, 2:end-1), [], 'all') - max(data(:, 2:end-1), [], 'all')*0.1;
    maxValue4 = min(data2(:, 2:end-1), [], 'all') + max(data(:, 2:end-1), [], 'all')*1.1;
    % Отображение данных перед началом цикла
    plotData(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4);
    
    % Получение гифки
    createGif(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4);
end

function [data, km] = loadData(filename)
    data = dlmread(filename, ';', 0, 0);
    km = 0:0.1:100;
end

function plotData(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4)
    figure;
    % Первый подграфик
    subplot(4, 1, 1);
    plot(km, data(1, 2:end-1), 'Color', 'b', LineWidth=2);
    hold on;
    xlabel('Труба, км');
    ylabel('Разница давлений, Па');
    title('Профиль разницы давлений');
    xlim([0, 100]);
    ylim([minValue, maxValue]);

    % Второй подграфик
    subplot(4, 1, 2);
    plot(km, data2(1, 2:end-1), 'Color', 'b');
    hold on;
    xlabel('Труба, км');
    ylabel('Давление, Па');
    title('Профиль давления');
    xlim([0, 100]);
    ylim([minValue2, maxValue2]);

    % Третий подграфик
    subplot(4, 1, 3);
    plot(km, data3(1, 2:end-1), 'Color', 'b', LineWidth=2);
    xlabel('Труба, км');
    ylabel('Плотность, кг/м3');
    title('Профиль плотности');
    xlim([0, 100]);
    ylim([minValue3, maxValue3]);
    
    % Четвертый подграфик
    subplot(4, 1, 4);
    plot(data2(1:end, 1), data2(1:end, end-1), 'Color', 'b', LineWidth=2);
    hold on;
    plot(data2(1,1), data2(1, end-1), "Marker",".","LineStyle","none",MarkerSize=20, Color='r');
    hold on;
    text(data2(1,1), data2(1, end-1), ['t = ' num2str(data2(1,1)) ', с'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off;
    xlabel('Время, с');
    ylabel('Давление, Па');
    title(['Времяной ряд давления на выходе']);
    xlim([0, 250000]);
    ylim([minValue4,maxValue4]);
    figure_size = [0, 0, 1920, 1080];
    set(gcf, 'Position', figure_size);
end

function createGif(data, data2, data3, km, minValue, maxValue, minValue2, maxValue2, minValue3, maxValue3, minValue4, maxValue4)
    % Получение кадра для первого кадра
    frame = getframe(gcf);
    [im, map] = rgb2ind(frame.cdata, 256, 'nodither');
    im(1, 1, 1, size(data, 1)) = 0;

    % Анимация
    for i = 2:size(data(:,1))
        % Отображение данных
        subplot(4, 1, 1);
        plot(km, data(i, 2:end-1), 'Color', 'r', LineWidth=2);
        hold on;
        plot(km, data(1, 2:end-1), 'Color', 'b', LineWidth=2);
        xlim([0, 100]);
        ylim([minValue, maxValue]);
        xlabel('Труба, км');
        ylabel('Разница давлений, Па');
        title('Профиль разницы давлений');
        hold off;

        subplot(4, 1, 2);
        plot(km, data2(i, 2:end-1), 'Color', 'r');
        hold on;
        plot(km, data2(1, 2:end-1), 'Color', 'b');
        xlim([0, 100]);
        ylim([minValue2, maxValue2]);
        xlabel('Труба, км');
        ylabel('Давление, Па');
        title('Профиль давления');
        hold off;

        subplot(4, 1, 3);
        plot(km, data3(i, 2:end-1), 'Color', 'b', LineWidth=2);
        xlabel('Труба, км');
        ylabel('Плотность, кг/м3');
        title('Профиль плотности');
        xlim([0, 100]);
        ylim([minValue3, maxValue3]);
        
        subplot(4, 1, 4);
        plot(data2(1:end, 1), data2(1:end, end-1), 'Color', 'b', LineWidth=2);
        hold on
        plot(data2(i,1), data2(i, end-1), "Marker",".","LineStyle","none",MarkerSize=20, Color='r')
        hold on;
        text(data2(i,1), data2(i, end-1), ['t = ' num2str(data2(i,1)) ', с'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        hold off;
        xlabel('Время, с');
        ylabel('Давление, Па');
        title(['Времяной ряд давления на выходе']);
        xlim([0, 250000]);
        ylim([minValue4,maxValue4]);
        figure_size = [0, 0, 1920, 1080];
        set(gcf, 'Position', figure_size);
        % Получение кадра
        frame = getframe(gcf);
        % Добавляем кадр к гифке
        im(:, :, 1, i) = rgb2ind(frame.cdata, map, 'nodither');
    end

    % Сохранение гифки в файл
    %filename = 'импульс_плотности.gif';
    filename = 'скачек_плотности.gif';
    imwrite(im, map, filename, 'DelayTime', 0.02, 'LoopCount', inf);
    disp(['Гифка сохранена в файл: ' filename]);
end