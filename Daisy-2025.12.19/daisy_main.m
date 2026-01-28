% daisy_main.m
% Daisy 有机质模型主程序
% 2025-12-17 lc

clear; clc; close all;

% 1. 设置模拟参数与时间范围 
Target_StartDate = datetime('2000-01-01 00:00:00'); 
Target_EndDate   = datetime('2001-01-01 00:00:00');

SimConfig.dt = 1.0;            % 时间步长 (小时)
SimConfig.n_layers = 5;        % 30 层土壤 (每层 10cm)
SimConfig.layer_depth = 10;    % cm
SimConfig.input_C_yearly = 500; % 每年输入的碳量 (g C/m2) 用于初始化平衡


% 2. 环境驱动数据（读取、筛选与对齐） 
% 设置文件名
weather_file = 'Inboundaryhours2000-2024.xlsx'; 
fprintf('正在读取气象数据: %s ...\n', weather_file);

% 自动检测导入选项，保留列名
opts = detectImportOptions(weather_file);
opts.VariableNamingRule = 'preserve'; 
RawTable = readtable(weather_file, opts);

% --- A. 解析时间列 (假设第1列是时间) ---
try
    Time_Col = RawTable{:, 1}; % 读取第一列
    
    if isdatetime(Time_Col)
        Source_Dates = Time_Col;
    elseif iscell(Time_Col) || isstring(Time_Col)
        % 尝试转换字符串格式
        Source_Dates = datetime(Time_Col, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    else
        error('第一列无法识别为时间格式。');
    end
catch ME
    error('时间列解析失败: %s', ME.message);
end

% --- B. 根据设定的起止时间筛选数据 ---
% 创建掩码
time_mask = (Source_Dates >= Target_StartDate) & (Source_Dates <= Target_EndDate);

% 提取对应的数据片段
Filtered_Dates = Source_Dates(time_mask);
Raw_Temp       = RawTable{time_mask, 2}; % 第2列: 温度
Raw_Rain       = RawTable{time_mask, 3}; % 第3列: 降雨 (mm)

if isempty(Filtered_Dates)
    error('筛选后的数据为空！请检查 Target_StartDate 和 Target_EndDate 是否在文件时间范围内。');
end

% --- C. 处理缺失值与单位转换 ---
% 处理 NaN
Raw_Temp(isnan(Raw_Temp)) = 15.0; % 缺测填 15度
Raw_Rain(isnan(Raw_Rain)) = 0.0;  % 缺测填 0雨

% 生成驱动向量
Temp_surface  = Raw_Temp;           % [℃]
Precip_hourly = Raw_Rain / 10.0;    % [mm] -> [cm]

% --- D. 自动更新模拟步数 ---
SimConfig.n_steps = length(Temp_surface);

fprintf('时间段截取成功:\n  开始: %s\n  结束: %s\n', string(Filtered_Dates(1)), string(Filtered_Dates(end)));
fprintf('模拟总时长: %d 小时 (%d 步)\n', hours(Target_EndDate - Target_StartDate), SimConfig.n_steps);

% 3. 初始化
fprintf('正在初始化分层土壤属性...\n');

% 生成分层土壤物理参数 (N_layers x 1 向量)
Params.Soil = setup_soil_profile(SimConfig);

% 初始化生化参数
[BiochemParams, State] = daisy_init(SimConfig);

% 合并参数到 Params 结构体
Params.SMB = BiochemParams.SMB; 
Params.SOM = BiochemParams.SOM;
Params.AOM = BiochemParams.AOM;
Params.DOM = BiochemParams.DOM;
Params.ClayModel = BiochemParams.ClayModel;
Params.Transport = BiochemParams.Transport;

% 确保 Env 驱动可以访问到 Clay
GlobalClayProfile = Params.Soil.Clay;

% 4. 结果存储预分配 (拆分为独立变量)  
n_layers = SimConfig.n_layers; 

% 结果变量初始化
Results.DateTimes = Filtered_Dates; 

% 每个变量都是 double 类型的二维矩阵 [时间 x 层]
Results.DateTimes = Filtered_Dates; % 保存时间轴用于绘图

% --- C 库 (g C/m2) ---
Results.SMB1_C = zeros(SimConfig.n_steps, n_layers); % SMB Slow
Results.SMB2_C = zeros(SimConfig.n_steps, n_layers); % SMB Fast

Results.SOM1_C = zeros(SimConfig.n_steps, n_layers); % SOM Slow
Results.SOM2_C = zeros(SimConfig.n_steps, n_layers); % SOM Fast
Results.SOM3_C = zeros(SimConfig.n_steps, n_layers); % SOM Inert

Results.AOM1_C = zeros(SimConfig.n_steps, n_layers); % AOM Slow
Results.AOM2_C = zeros(SimConfig.n_steps, n_layers); % AOM Fast

% --- N 库 (g N/m2) ---
Results.SMB1_N = zeros(SimConfig.n_steps, n_layers);
Results.SMB2_N = zeros(SimConfig.n_steps, n_layers);

Results.SOM1_N = zeros(SimConfig.n_steps, n_layers);
Results.SOM2_N = zeros(SimConfig.n_steps, n_layers);
Results.SOM3_N = zeros(SimConfig.n_steps, n_layers);

Results.AOM1_N = zeros(SimConfig.n_steps, n_layers);
Results.AOM2_N = zeros(SimConfig.n_steps, n_layers);

% --- 矿质氮 ---
Results.Soil_NH4 = zeros(SimConfig.n_steps, n_layers);
Results.Soil_NO3 = zeros(SimConfig.n_steps, n_layers);
Results.Leaching_NO3 = zeros(SimConfig.n_steps, 1); 
Results.Leaching_NH4 = zeros(SimConfig.n_steps, 1); % 淋溶流失

% 累积呼吸
Results.CO2_Cum = zeros(SimConfig.n_steps, 1);
Current_CO2_Cum = 0;

% 5. 生成土壤水热剖面
% 构造输入结构体
WeatherInput.Precip_hourly = Precip_hourly;
WeatherInput.Temp_surface  = Temp_surface;

% 确保 daisy_drive_soil_richards.m 在路径中
EnvDrivers = daisy_drive_soil_richards(WeatherInput, SimConfig, Params.Soil);

% 提取结果
Temp_profile     = EnvDrivers.Temp_profile;
Theta_profile    = EnvDrivers.Theta_profile;
Moisture_profile = EnvDrivers.Moisture_profile;
Water_Flux       = EnvDrivers.Water_Flux;
Clay_profile     = EnvDrivers.Clay;

% 6. 主循环 
fprintf('开始模拟 (%d 步)...\n', SimConfig.n_steps);
tic;

for t = 1:SimConfig.n_steps
    % A. 准备环境
    Env.Temp = Temp_profile(t, :)';
    Env.h = Moisture_profile(t, :)';
    Env.Clay = Clay_profile;
    Env.Theta_new = Theta_profile(t, :)';

    if t == 1
        Env.Theta_old = Theta_profile(t, :)';
    else
        Env.Theta_old = Theta_profile(t-1, :)';
    end
    Env.q_edge = Water_Flux(t, :)'; % 界面通量 [N+1 x 1]

    % B. 执行步长
    [State, Fluxes] = daisy_step(State, Params, Env, SimConfig.dt);
    
    % C. 详细记录结果 (分库记录)
    
    % State.SMB.C 是 [层数 x 库数]，转置或按列取使其适配 [1 x 层数]
    % 注意：State.SMB.C(:, 1) 是第一列(SMB1)的所有层数据
    
    % --- 保存 SMB ---
    Results.SMB1_C(t, :) = State.SMB.C(:, 1)'; 
    Results.SMB2_C(t, :) = State.SMB.C(:, 2)';
    Results.SMB1_N(t, :) = State.SMB.N(:, 1)'; 
    Results.SMB2_N(t, :) = State.SMB.N(:, 2)';
    
    % --- 保存 SOM ---
    Results.SOM1_C(t, :) = State.SOM.C(:, 1)'; 
    Results.SOM2_C(t, :) = State.SOM.C(:, 2)'; 
    Results.SOM3_C(t, :) = State.SOM.C(:, 3)';
    Results.SOM1_N(t, :) = State.SOM.N(:, 1)'; 
    Results.SOM2_N(t, :) = State.SOM.N(:, 2)'; 
    Results.SOM3_N(t, :) = State.SOM.N(:, 3)';
    
    % --- 保存 AOM ---
    Results.AOM1_C(t, :) = State.AOM.C(:, 1)'; 
    Results.AOM2_C(t, :) = State.AOM.C(:, 2)';
    Results.AOM1_N(t, :) = State.AOM.N(:, 1)'; 
    Results.AOM2_N(t, :) = State.AOM.N(:, 2)';
    
    % --- 保存 矿质氮 ---
    Results.Soil_NH4(t, :) = State.Soil.NH4';
    Results.Soil_NO3(t, :) = State.Soil.NO3';

    if isfield(Fluxes, 'Leaching')
        Results.Leaching_NO3(t) = Fluxes.Leaching.NO3;
        Results.Leaching_NH4(t) = Fluxes.Leaching.NH4;
    end
    
    % 累积呼吸
    Current_CO2_Cum = Current_CO2_Cum + sum(Fluxes.CO2_Total);
    Results.CO2_Cum(t) = Current_CO2_Cum;
    
   % D. 简单的农事操作示例 (按日期触发)
    % 检查当前日期
    CurrentDate = Results.DateTimes(t);

    % 示例：每年 10月1日 添加秸秆
    if month(CurrentDate) == 10 && day(CurrentDate) == 1 && hour(CurrentDate) == 0
        fprintf('[%s] 执行秸秆还田...\n', string(CurrentDate));
        added_C = 200; % g/m2
        added_N = added_C / 20.0;
        State.AOM.C(1, 2) = State.AOM.C(1, 2) + added_C; % 加到表层 AOM2
        State.AOM.N(1, 2) = State.AOM.N(1, 2) + added_N;
    end

    % 示例：每年 10月15日 耕作 (混合表层)
    if month(CurrentDate) == 10 && day(CurrentDate) == 15 && hour(CurrentDate) == 0
        fprintf('[%s] 执行耕作混合...\n', string(CurrentDate));
        mix = 1:2; % 混合前两层
        State.SOM.C(mix,:) = repmat(mean(State.SOM.C(mix,:),1), 2, 1);
        % ... (其他库也混合，此处略)
    end
end
toc;

% 7. 绘图
n_layers = size(Results.SMB1_C, 2);
layer_colors = parula(n_layers); 

figure('Name', 'Detailed Carbon Pools Dynamics', 'Color', 'w', 'Position', [100, 100, 1600, 1000]);

pool_vars = {'SMB1_C', 'SMB2_C', 'AOM1_C', 'AOM2_C', 'SOM1_C', 'SOM2_C', 'SOM3_C'};
pool_titles = {
    'SMB1 (Slow Microbes)', 'SMB2 (Fast Microbes)', ...
    'AOM1 (Slow Added Matter)', 'AOM2 (Fast Added Matter)', ...
    'SOM1 (Slow Humus)', 'SOM2 (Fast Humus)', 'SOM3 (Inert Humus)'
};

for i = 1:7
    subplot(3, 3, i); 
    hold on;
    data = Results.(pool_vars{i});
    
    for z = 1:n_layers
        % 使用 Results.DateTimes 作为 X 轴
        plot(Results.DateTimes, data(:, z), ... 
             'Color', layer_colors(z, :), ...
             'LineWidth', 1.5);
    end
    
    title(pool_titles{i}, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('g C / m^2');
    grid on;
    % 自动调整日期格式
    xtickformat('yyyy-MM'); 
    
    if i == 1
        legend(strcat('Layer ', string(1:n_layers)), 'Location', 'best');
    end
end

% 绘制温度和水分辅助图
subplot(3, 3, 8);
plot(Results.DateTimes, Temp_surface, 'r-');
title('Surface Temperature'); ylabel('°C');
xtickformat('yyyy-MM'); grid on;

subplot(3, 3, 9);
plot(Results.DateTimes, Precip_hourly * 10, 'b-'); % cm -> mm
title('Rainfall'); ylabel('mm/h');
xtickformat('yyyy-MM'); grid on;

sgtitle(['Simulation: ' char(string(Target_StartDate)) ' to ' char(string(Target_EndDate))]);