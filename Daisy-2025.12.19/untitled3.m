% ========================================================
% 多层土壤诊断脚本 (5 Layers Diagnosis)
% 适用：检查各层养分、水分及有机质分解差异
% ========================================================

figure('Name', 'Multi-Layer Diagnosis', 'Color', 'w', 'Position', [50, 50, 1400, 900]);

% --- 绘图通用设置 ---
layers_to_plot = 1:5;             % 要画的层数
colors = lines(5);                % 生成5种对比鲜明的颜色
legend_str = cell(1, 5);
for i = 1:5, legend_str{i} = ['Layer ' num2str(i)]; end

% ========================================================
% 1. 铵态氮 (NH4) - 5层对比
% ========================================================
subplot(2,2,1);
hold on;
if isfield(Results, 'Soil_NH4')
    for i = layers_to_plot
        if size(Results.Soil_NH4, 2) >= i
            plot(Results.DateTimes, Results.Soil_NH4(:, i), ...
                 'LineWidth', 1.5, 'Color', colors(i,:));
        end
    end
    title('土壤铵态氮 (NH_4) 动态');
    ylabel('g N / m^2');
    legend(legend_str, 'Location', 'Best');
else
    text(0.5,0.5,'未找到 NH4 数据','HorizontalAlignment','center');
end
grid on; hold off;

% ========================================================
% 2. 硝态氮 (NO3) - 5层对比
% ========================================================
subplot(2,2,2);
hold on;
if isfield(Results, 'Soil_NO3')
    for i = layers_to_plot
        if size(Results.Soil_NO3, 2) >= i
            plot(Results.DateTimes, Results.Soil_NO3(:, i), ...
                 'LineWidth', 1.5, 'Color', colors(i,:));
        end
    end
    title('土壤硝态氮 (NO_3) 动态');
    ylabel('g N / m^2');
    legend(legend_str, 'Location', 'Best');
else
    text(0.5,0.5,'未找到 NO3 数据','HorizontalAlignment','center');
end
grid on; hold off;

% ========================================================
% 3. 土壤水分 (Theta) - 5层对比
% ========================================================
subplot(2,2,3);
hold on;
% 尝试查找常见的水分变量名
Data_Theta = [];
if exist('Theta_profile', 'var')
    Data_Theta = Theta_profile;
elseif isfield(Results, 'Theta')
    Data_Theta = Results.Theta;
elseif isfield(Results, 'Soil_Water')
    Data_Theta = Results.Soil_Water;
end

if ~isempty(Data_Theta)
    for i = layers_to_plot
        if size(Data_Theta, 2) >= i
            plot(Results.DateTimes, Data_Theta(:, i), ...
                 'LineWidth', 1.5, 'Color', colors(i,:));
        end
    end
    title('土壤体积含水量 (\theta)');
    ylabel('cm^3 / cm^3');
    ylim([0 0.6]); % 通常土壤水分在 0-0.5 之间，可根据需要调整
    legend(legend_str, 'Location', 'Best');
else
    text(0.5,0.5,'未找到水分变量 (Theta_profile)','HorizontalAlignment','center');
end
grid on; hold off;

% ========================================================
% 4. 添加有机质 (AOM2) - 5层对比
% ========================================================
subplot(2,2,4);
hold on;
% 这里默认画 AOM2 (Fast Added Matter)，因为它是最活跃的
% 如果您想看总和，可以改为 Results.AOM1_C + Results.AOM2_C
if isfield(Results, 'AOM2_C')
    for i = layers_to_plot
        if size(Results.AOM2_C, 2) >= i
            plot(Results.DateTimes, Results.AOM2_C(:, i), ...
                 'LineWidth', 1.5, 'Color', colors(i,:));
        end
    end
    title('外源有机质残留 (AOM2)');
    ylabel('g C / m^2');
    legend(legend_str, 'Location', 'Best');
else
    text(0.5,0.5,'未找到 AOM2 数据','HorizontalAlignment','center');
end
grid on; hold off;