function SoilParams = setup_soil_profile(SimConfig)
% SETUP_SOIL_PROFILE 定义分层土壤属性并映射到计算网格
%
% 依据: Daisy 通常根据土壤图 (Soil Map) 定义 Horizon (A层, B层, C层...)
%       这里我们手动定义一个典型的农业土壤剖面 (例如 Sandy Loam)
% 输入: SimConfig 结构体, 包含模拟层数 n_layers 和层厚 layer_depth 
% 输出: SoilParams 结构体，包含每一层网格的物理参数向量 (N_layers x 1)

% 2025-12-16 lc

% 从配置中提取层数和层厚
    n_layers = SimConfig.n_layers;
    layer_depth = SimConfig.layer_depth; % cm
    
    %  1. 定义发生层 (Horizons) 
    % 这里定义 3 个典型的土层
    
    % Horizon 1: 表层/耕作层 (0 - 30 cm) - 有机质较高，疏松
    H(1).BottomDepth = 30;
    H(1).Clay   = 0.12; % 粘土含量 12%
    H(1).Sand   = 0.60;
    H(1).OM     = 0.025; % 有机质 2.5%
    H(1).BD     = 1.35;  % 容重 (g/cm3)
    % 水力参数 (van Genuchten - Mualem)
    H(1).alpha  = 0.024; % 1/cm
    H(1).n      = 1.42;
    H(1).K_sat  = 12.0;  % cm/h (表层导水率高)
    H(1).Th_sat = 0.43;  % 饱和含水量
    H(1).Th_res = 0.01;  % 残余含水量

    % Horizon 2: 犁底层/亚表层 (30 - 60 cm) - 压实，导水率低
    H(2).BottomDepth = 60;
    H(2).Clay   = 0.15; 
    H(2).Sand   = 0.55;
    H(2).OM     = 0.01;
    H(2).BD     = 1.55;  % 容重较大
    % 水力参数
    H(2).alpha  = 0.018;
    H(2).n      = 1.35;
    H(2).K_sat  = 2.5;   % 导水率显著降低
    H(2).Th_sat = 0.40;
    H(2).Th_res = 0.01;

    % Horizon 3: 母质层/深层 (60 - 底部)
    H(3).BottomDepth = 9999;
    H(3).Clay   = 0.18;
    H(3).Sand   = 0.50;
    H(3).OM     = 0.002;
    H(3).BD     = 1.60;
    % 水力参数
    H(3).alpha  = 0.015;
    H(3).n      = 1.30;
    H(3).K_sat  = 1.0;   % 粘土增多，导水慢
    H(3).Th_sat = 0.38;
    H(3).Th_res = 0.01;

    % --- 2. 映射到计算网格 (Discretization) ---
    % 初始化向量
    SoilParams.Clay   = zeros(n_layers, 1);
    SoilParams.Sand   = zeros(n_layers, 1);
    SoilParams.OM     = zeros(n_layers, 1);
    SoilParams.BD     = zeros(n_layers, 1);
    SoilParams.alpha  = zeros(n_layers, 1);
    SoilParams.n      = zeros(n_layers, 1);
    SoilParams.m      = zeros(n_layers, 1); % m = 1 - 1/n
    SoilParams.K_sat  = zeros(n_layers, 1);
    SoilParams.Th_sat = zeros(n_layers, 1);
    SoilParams.Th_res = zeros(n_layers, 1);
    
    for i = 1:n_layers
        % 计算当前层中心深度
        z_mid = (i - 5) * layer_depth;
        
        % 查找所属的 Horizon
        found = false;
        for k = 1:length(H)
            if z_mid <= H(k).BottomDepth
                sel = k;
                found = true;
                break;
            end
        end
        if ~found, sel = length(H); end % 超过最深层使用最后一层属性
        
        % 赋值
        SoilParams.Clay(i)   = H(sel).Clay;
        SoilParams.Sand(i)   = H(sel).Sand;
        SoilParams.OM(i)     = H(sel).OM;
        SoilParams.BD(i)     = H(sel).BD;
        SoilParams.alpha(i)  = H(sel).alpha;
        SoilParams.n(i)      = H(sel).n;
        SoilParams.m(i)      = 1 - 1/H(sel).n;
        SoilParams.K_sat(i)  = H(sel).K_sat;
        SoilParams.Th_sat(i) = H(sel).Th_sat;
        SoilParams.Th_res(i) = H(sel).Th_res;
    end
    
    fprintf('土壤剖面构建完成: %d 层, 分为 %d 个发生层。\n', n_layers, length(H));
end