function State = daisy_init_solver(Params, Config)
% 使用线性方程组求解初始稳态 (Gaussian Elimination)
% 对应 src/daisy/organic_matter/organic_std.C::partition
% 2025-12-17 lc
    
    n_lay = Config.n_layers;
    
    % 1. 确定平均输入速率 (g C/m2/h)
    Input_Rate = Config.input_C_yearly / (365 * 24);
    
    % 2. 构建系数矩阵 M * X = B
    % 变量 X = [SMB1, SMB2, SOM1, SOM2]' (SOM3 是惰性的，手动设定)
    M = zeros(4, 4);
    B = zeros(4, 1);
    
    % 假设输入流向 SMB2 (来自 AOM 分解)
    % AOM -> SMB2 = Input * Eff_AOM * Frac_AOM_SMB2
    B(2) = Input_Rate * 0.5 * 1.0; 
    
    % 填充矩阵 (简化版: 忽略维持呼吸对稳态的非线性影响，仅考虑 Turnover)
    % 这里我们构建 C 的流向平衡方程: In - Out = 0
    
    % SMB1 方程 (Index 1)
    M(1,1) = -Params.SMB.k_turnover(1); % Out
    M(1,2) = Params.SMB.k_turnover(2) * Params.SMB.efficiency(2) * Params.SMB.fractions(2,1); % In from SMB2
    M(1,3) = Params.SOM.k_turnover(1) * Params.SOM.efficiency(1) * Params.SOM.fractions(1,1); % In from SOM1
    M(1,4) = Params.SOM.k_turnover(2) * Params.SOM.efficiency(2) * Params.SOM.fractions(2,1); % In from SOM2
    
    % SMB2 方程 (Index 2)
    M(2,2) = -Params.SMB.k_turnover(2);
    M(2,1) = Params.SMB.k_turnover(1) * Params.SMB.efficiency(1) * Params.SMB.fractions(1,2);
    M(2,3) = Params.SOM.k_turnover(1) * Params.SOM.efficiency(1) * Params.SOM.fractions(1,2);
    M(2,4) = Params.SOM.k_turnover(2) * Params.SOM.efficiency(2) * Params.SOM.fractions(2,2);
    
    % SOM1 方程 (Index 3)
    M(3,3) = -Params.SOM.k_turnover(1);
    M(3,1) = Params.SMB.k_turnover(1) * Params.SMB.efficiency(1) * Params.SMB.fractions(1,3);
    M(3,2) = Params.SMB.k_turnover(2) * Params.SMB.efficiency(2) * Params.SMB.fractions(2,3);
    M(3,4) = Params.SOM.k_turnover(2) * Params.SOM.efficiency(2) * Params.SOM.fractions(2,3);
    
    % SOM2 方程 (Index 4)
    M(4,4) = -Params.SOM.k_turnover(2);
    M(4,1) = Params.SMB.k_turnover(1) * Params.SMB.efficiency(1) * Params.SMB.fractions(1,4);
    M(4,2) = Params.SMB.k_turnover(2) * Params.SMB.efficiency(2) * Params.SMB.fractions(2,4);
    M(4,3) = Params.SOM.k_turnover(1) * Params.SOM.efficiency(1) * Params.SOM.fractions(1,4);
    
    % 求解 X = M \ -B (因为 M*X + B = 0 -> M*X = -B)
    % X 的单位是 g C / m2
    X = M \ (-B);
    
    % 3. 分配到各层 (按深度指数衰减)
    State.SMB.C = zeros(n_lay, 2); State.SOM.C = zeros(n_lay, 3);
    
    for z = 1:n_lay
        depth_factor = exp(-(z-1)*0.5); % 模拟分布
        
        State.SMB.C(z, 1) = X(1) * depth_factor ; 
        State.SMB.C(z, 2) = X(2) * depth_factor ;
        State.SOM.C(z, 1) = X(3) * depth_factor ;
        State.SOM.C(z, 2) = X(4) * depth_factor ;
        State.SOM.C(z, 3) = 2000 * depth_factor ; % Inert 手动指定（g/m2）
    end
    
    % 4. 初始化 N 和其他
    State.SMB.N = State.SMB.C ./ Params.SMB.CN;
    State.SOM.N = State.SOM.C ./ Params.SOM.CN;
    
    State.AOM.C = zeros(n_lay, 2); State.AOM.N = zeros(n_lay, 2);
    State.DOM.C = zeros(n_lay, 1); State.DOM.N = zeros(n_lay, 1);
    State.DOM.C_Sorb = zeros(n_lay, 1); State.DOM.N_Sorb = zeros(n_lay, 1); % 吸附态 DOM
    
    State.Soil.NH4 = ones(n_lay, 1) * 5; % g N/m2
    State.Soil.NO3 = ones(n_lay, 1) * 10;
    
    disp('初始化求解完成。');
end