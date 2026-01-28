function [Params, State] = daisy_init(Config)
% 初始化 Daisy 模型的参数，并调用求解器计算初始稳态
% 输入：Config（配置结构体，通常包含模拟年份、土层厚度等基本设置）
% 输出：Params（包含所有物理化学参数的结构体）
%       State（包含所有状态变量初始值的结构体）

% 2025-12-16 lc
   
    % 参数定义 (Params) 
    % 速率单位: 1/h
    
    % 1. SMB (Microbial Biomass) - 2个库: [Slow, Fast]
    Params.SMB.k_turnover = [7.7e-6, 4.2e-4]; % k_turnover: 微生物的死亡/周转速率
    Params.SMB.k_maint    = [7.5e-5, 4.2e-4]; % k_maint: 微生物的维持呼吸速率，这部分碳直接变成CO2
    Params.SMB.efficiency = [0.6, 0.6]; % efficiency: 碳利用效率
    Params.SMB.CN         = [6.7, 6.7]; % CN: 微生物碳氮比

    % 1=SMB1（慢菌）
    % 2=SMB2（快菌）
    % 3=SOM1（Slow 较难分解的结构性有机质）
    % 4=SOM2（Fast 较易分解的活性腐殖质）
    % 5=SOM3（Inert 惰性碳库，模型时间尺度内不分解） 
    % 分配矩阵 Fractions: [Source x Target]
    % Daisy模型的核心逻辑之一，模拟微生物残体向土壤有机质库的转化
    % 第1行：SMB1(慢菌)死亡，60% 变成 SMB2，40% 变成 SOM2（腐殖质）
    % 第2行：SMB2(快菌)死亡，40% 变成 SMB2，60% 变成 SOM2（腐殖质）
    Params.SMB.fractions = [
        0.0, 0.6, 0.0, 0.4, 0.0;
        0.0, 0.4, 0.0, 0.6, 0.0
    ];

    % 2. SOM (Soil Organic Matter) - 3个库: [Slow, Fast, Inert]
    Params.SOM.k_turnover = [1.8e-6, 5.8e-6, 0.0];
    Params.SOM.efficiency = [0.4, 0.5, 0.0];
    Params.SOM.CN         = [10.0, 10.0, 10.0];
    %第1行: SOM1 分解后，100% 流向 SMB1（成为慢速微生物的食物）。
    %第2行: SOM2 分解后，70% 流向 SMB1，30% 转变为 SOM1（化学固化/聚合作用）。
    %第3行: 惰性库（虽然不分解，但保留格式一致性）。
    Params.SOM.fractions = [
        1.0, 0.0, 0.0, 0.0, 0.0;
        0.7, 0.0, 0.3, 0.0, 0.0;
        0.0, 0.0, 0.0, 0.0, 1.0
    ];

    % 3. AOM (Added Organic Matter) - 2个库: [Slow, Fast]
    % 模拟施肥或者秸秆还田
    Params.AOM.k_turnover = [2.0e-4, 2.0e-3];
    Params.AOM.efficiency = [0.5, 0.5];
    Params.AOM.CN         = [50.0, 20.0]; % 高C/N，代表新鲜植物残体，导致微生物分解时需要从土壤中固氮
    % AOM 产物均流向 SMB2，模拟新鲜有机物引发的微生物爆发（激发效应）
    Params.AOM.fractions = [
        0.0, 1.0, 0.0, 0.0, 0.0; 
        0.0, 1.0, 0.0, 0.0, 0.0
    ];
    
    % 4. DOM (Dissolved Organic Matter)
    Params.DOM.k_turnover = 0.05;  % 周转率 (1/h)，通常比 SOM 快
    Params.DOM.efficiency = 0.6;   % 微生物利用效率
    Params.DOM.fractions  = [0.0, 1.0]; % DOM 分解后 100% 进入 SMB2 (Fast)

    % 5. 吸附模型配置
    % 定义吸附的数学形式
    % 'Linear'，吸附量与浓度成正比（简单，但在高浓度下不真实）。
    % 'Langmuir'（当前设定），假设土壤吸附位点有限，浓度高时吸附会饱和（更符合物理真实）
    Params.DOM.SorptionModel = 'Langmuir'; % 可选: 'Linear', 'Langmuir', 'Freundlich'
    Params.DOM.alpha = 0.1;        % 吸附动力学速率
    % Langmuir 参数
    Params.DOM.K_adsorption = 0.1; % Langmuir 方程中的亲和力常数（K_L）
    Params.DOM.S_max = 5e-4; % 最大吸附容量    
    % Linear 参数 (如果选 Linear)
    Params.DOM.K_dist = 5.0; % 线性分配系数(K_d)

    % 6. porosity 孔隙度
    Params.Transport.porosity = 0.5; %通常约 0.4-0.5

    % 7. 溶质运移参数
    Params.Transport.dispersivity = 5.0; % [cm] 水动力弥散度 (λ) 值越大溶质锋面越模糊（平滑）
    Params.Transport.diffusion_coeff = 1e-5 * 3600; % [cm2/h] 分子扩散系数（D0）
    Params.Transport.C_rain = 0.0; % [g/cm3] 边界条件，降雨中的溶质浓度，这里假设降雨不含DOM或氮素
    
    % 8. 粘土效应选择 ---
    % 'Standard': 经典 Daisy (粘土保护SOM, 增加SMB维持消耗)
    % 'Biomod'  : BIOMOD 模型 (粘土抑制SMB分解和维持, 对SOM无影响)
    % 'None'    : 不考虑粘土效应
    Params.ClayModel = 'Standard';

    

    % 求解初始碳氮库大小（调用 daisy_int_solver）
    State = daisy_init_solver(Params, Config);
end