function EnvDrivers = daisy_drive_soil_richards(Weather, SimConfig, SoilParams)
% DAISY_DRIVE_SOIL_RICHARDS (Layered Version)
% 基于分层 Richards 方程和动态热参数的土壤物理驱动
%
% 输入:
%   Weather:    包含 Precip_hourly, Temp_surface
%   SimConfig:  配置
%   SoilParams: 由 setup_soil_profile 生成的结构体 (包含各层的向量参数)
% 输出：
%   EnvDrivers（T, Theta_new, q） 

% 2025-12-16 lc

    fprintf('正在计算分层土壤物理过程 (Richards + Dynamic Heat)...\n');

% 提取配置
    N_steps = SimConfig.n_steps; % 总时长
    N_layers = SimConfig.n_layers; % 土层
    dz = SimConfig.layer_depth; % 层厚
    
    % 1. 建立全 0 矩阵
    Temp_prof  = zeros(N_steps, N_layers);
    Theta_prof = zeros(N_steps, N_layers);
    h_prof     = zeros(N_steps, N_layers);
    Flux_prof  = zeros(N_steps, N_layers + 1); % 界面比土层多1
    
    % 2. 初始化状态
    % 设定初始水势: -100 cm (田间持水量附近，湿润但不饱和)
    h_current = ones(N_layers, 1) * -100; 
    % 假设整个土层的初始温度都等于第一个小时的地表气温（简化假设）
    T_current = ones(N_layers, 1) * Weather.Temp_surface(1);
    
    % 提取土壤属性向量 (简化结构体字段)
    S_alpha = SoilParams.alpha;
    S_n     = SoilParams.n;
    S_m     = SoilParams.m;
    S_Ksat  = SoilParams.K_sat;
    S_Th_s  = SoilParams.Th_sat;
    S_Th_r  = SoilParams.Th_res;
    
    % 计算土壤固体部分的热容
    % 简化计算固体热容 (J/cm3/K): 矿物~2.0, 有机质~2.5
    C_solid_vec = 2.0 * (1 - SoilParams.OM) + 2.5 * SoilParams.OM;
    
    % 3. 时间步循环
    for t = 1:N_steps
        
        q_rain = Weather.Precip_hourly(t); 
        T_surf = Weather.Temp_surface(t);
        
        % 如果时间在早上6点到晚上18点之间，且没有下雨，则设定一个固定的蒸发速率 0.02 cm/h。
        % 否则为 0。这是一个简化的物理过程
        ET_pot = 0.0;
        hour_of_day = mod(t, 24);
        if hour_of_day > 6 && hour_of_day < 18 && q_rain == 0
            ET_pot = 0.02; % cm/h
        end
        
        dt_hour = 1.0; 
        time_done = 0;
        
        % 自适应子步长循环
        while time_done < dt_hour
            % 根据降雨强度调整步长
            if q_rain > 1.0, dt = 0.05; 
            elseif q_rain > 0, dt = 0.1;
            else, dt = 0.25; 
            end
            
            if time_done + dt > dt_hour, dt = dt_hour - time_done; end
            

            % A. 求解 Richards 方程（水分运动）
            h_old_sub = h_current;
            h_iter = h_current;
            
            % Picard 迭代循环 (处理非线性)
            for iter = 1:15
                % 1. 计算每一层的导水率 K 和容水度 C (传入向量参数)
                [K_vec, C_vec] = calc_hydraulic_layered(h_iter, S_alpha, S_n, S_m, S_Ksat, S_Th_s, S_Th_r);
                
                % 2. 组装矩阵，初始化三对角矩阵的系数向量和右端项向量，用于解线性方程组
                A = zeros(N_layers, 1); B = zeros(N_layers, 1);
                C_mtx = zeros(N_layers, 1); D = zeros(N_layers, 1);
                
                % 建立差分方程
                for i = 1:N_layers
                    % 几何平均导水率
                    % 计算当前层与上层（K_up）及下层（K_dn）界面处的导水率
                    % 采用几何平均算法（sqrt(K1*K2)），这在处理层间差异大时更稳定
                    if i > 1, K_up = sqrt(K_vec(i)*K_vec(i-1)); else, K_up = K_vec(i); end
                    if i < N_layers, K_dn = sqrt(K_vec(i)*K_vec(i+1)); else, K_dn = K_vec(i); end
                    
                    % 计算地表净流量：降雨-蒸发
                    q_surf_net = q_rain - ET_pot;
                    
                    % 特殊处理上边界条件，将净流量 q_surf_net 加入方程的右端项 D
                    if i == 1
                        % 上边界: Flux = q_surf
                        A(i) = 0;
                        B(i) = C_vec(i)/dt + K_dn/dz;
                        C_mtx(i) = -K_dn/dz; % 注意这里移项后的符号
                        D(i) = C_vec(i)/dt * h_old_sub(i) + q_surf_net/dz - K_dn/dz; 
                    % 特殊处理下边界条件，假设为自由排水（重力流），即没有吸力梯度，水自然流出
                    elseif i == N_layers
                        % 下边界: 自由排水 Flux_out = K_i
                        A(i) = -K_up/dz;
                        B(i) = C_vec(i)/dt + K_up/dz;
                        C_mtx(i) = 0;
                        D(i) = C_vec(i)/dt * h_old_sub(i) + K_up/dz - K_vec(i)/dz;
                    % 标准的隐式差分格式。A关联上层，C_mtx关联下层，B是本层系数
                    else
                        % 内部节点
                        A(i) = -K_up/dz;
                        B(i) = C_vec(i)/dt + (K_up + K_dn)/dz;
                        C_mtx(i) = -K_dn/dz;
                        D(i) = C_vec(i)/dt * h_old_sub(i) + (K_up - K_dn)/dz;
                    end
                end
                
                % 调用追赶法求解器tridiag_solver，解出线性方程组，得到这一步新的水势 h_new_sub
                h_new_sub = tridiag_solver(A, B, C_mtx, D);
                
                % 收敛检查
                if max(abs(h_new_sub - h_iter)) < 1e-2
                    break;
                end
                h_iter = h_new_sub;
            end
            h_current = h_iter;
            
            % B. 求解热传输 (动态热参数)
            % 1. 根据刚才算出的新水势 h_current，反推当前的准确含水量 Theta_curr 和导水率 K_curr (用于热参数)
            [Theta_curr, K_curr, ] = calc_hydraulic_layered(h_current, S_alpha, S_n, S_m, S_Ksat, S_Th_s, S_Th_r);
            
            % 2. 计算每一层的热参数 (Daisy Physics)
            % 计算土壤体积热容 固体热容×固体体积比例+水的热容（4.18）×含水量
            % C_soil = f_solid*C_solid + theta*C_water
            % C_water ~ 4.18 J/cm3/K
            C_soil_vol = C_solid_vec .* (1 - SoilParams.Th_sat) + Theta_curr * 4.18;
            
            % 计算土壤导热率 Lambda (Kersten Number method approx)
            % 简化: Lambda = A + B * Theta (或者更复杂的)
            % 这里使用简单的 Campbell 经验公式模拟沙壤土
            Lambda_vec = 15 + 150 * (Theta_curr ./ SoilParams.Th_sat); % [J/cm/h/K] 估算值
            
            % 3. 利用达西定律计算各层界面上的水流速度 
            % 为了计算热对流对流热（水流动带走的热量）
            %q_edges(1) 是雨水入渗。
            % dh_dz - 1 是总水势梯度（含重力项 -1）
            q_edges = zeros(N_layers+1, 1);
            q_edges(1) = q_rain;
            for i = 1:N_layers-1
                K_avg = sqrt(K_curr(i)*K_curr(i+1));
                dh_dz = (h_current(i+1)-h_current(i))/dz;
                q_edges(i+1) = -K_avg * (dh_dz - 1); 
            end
            q_edges(end) = K_curr(end);
            
            % 4. 求解热方程
            % 调用 solve_heat_layered 函数解热传导+对流方程，算出新温度
            T_new_sub = solve_heat_layered(T_current, T_surf, q_edges, Lambda_vec, C_soil_vol, dt, dz, N_layers);
            T_current = T_new_sub;
            
            time_done = time_done + dt;
        end
        
        % 存储结果
        h_prof(t, :)     = h_current;
        Theta_prof(t, :) = calc_theta_layered(h_current, S_alpha, S_n, S_m, S_Th_s, S_Th_r);
        Temp_prof(t, :)  = T_current;
        Flux_prof(t, :)  = q_edges;
    end
    
    EnvDrivers.Temp_profile     = Temp_prof;
    EnvDrivers.Theta_profile    = Theta_prof;
    EnvDrivers.Moisture_profile = h_prof;
    EnvDrivers.Water_Flux       = Flux_prof;
    EnvDrivers.Clay             = SoilParams.Clay; % 传递分层粘土数据
    
    fprintf('土壤物理过程计算完成。\n');
end

% 辅助函数: 向量化 van Genuchten
% 计算 van Genuchten 水分特征曲线
function [theta, K, C] = calc_hydraulic_layered(h, alpha, n, m, K_sat, Th_s, Th_r)
    h = min(h, -1e-3);
    rel_term = (1 + abs(alpha .* h).^n).^m;
    theta = Th_r + (Th_s - Th_r) ./ rel_term;
    
    % 计算容水度C (即 dTheta/dh)，这是 Richards 方程求解所需的微分项
    num = alpha .* n .* m .* (Th_s - Th_r) .* (abs(alpha .* h).^(n-1));
    den = (1 + abs(alpha .* h).^n).^(m+1);
    C = num ./ den;
    
    % 计算 Mualem-van Genuchten 导水率 K。这决定了水在土壤里流动的快慢
    Se = (theta - Th_r) ./ (Th_s - Th_r);
    Se = max(0, min(1, Se));
    term_k = (1 - (1 - Se.^(1./m)).^m).^2;
    K = K_sat .* (Se.^0.5) .* term_k;
end

% 简化版的辅助函数，只返回含水量，用于最后结果输出
function theta = calc_theta_layered(h, alpha, n, m, Th_s, Th_r)
    h = min(h, -1e-3);
    rel_term = (1 + abs(alpha .* h).^n).^m;
    theta = Th_r + (Th_s - Th_r) ./ rel_term;
end

% --- 辅助函数: 分层热方程求解 ---
function T_new = solve_heat_layered(T_old, T_surf, q_edge, Lambda, C_vol, dt, dz, N)
    Cw = 4.18; % 水的热容
    dT_dt = zeros(N, 1);
    
    for i = 1:N
        if i == 1, T_up = T_surf; else, T_up = T_old(i-1); end
        if i == N, T_dn = T_old(i); else, T_dn = T_old(i+1); end
        
        % 导热率平均
        if i==1, L_up=Lambda(i); else, L_up=(Lambda(i)+Lambda(i-1))/2; end
        if i==N, L_dn=Lambda(i); else, L_dn=(Lambda(i)+Lambda(i+1))/2; end
        
        % 传导通量
        J_cond_in  = L_up * (T_up - T_old(i)) / dz;
        J_cond_out = L_dn * (T_old(i) - T_dn) / dz;
        
        % 对流通量
        q_in = q_edge(i); q_out = q_edge(i+1);
        T_flow_in  = (q_in > 0) * T_up + (q_in <= 0) * T_old(i);
        T_flow_out = (q_out > 0) * T_old(i) + (q_out <= 0) * T_dn;
        J_conv_in  = q_in * Cw * T_flow_in;
        J_conv_out = q_out * Cw * T_flow_out;
        
        dT_dt(i) = (J_cond_in - J_cond_out + J_conv_in - J_conv_out) / (C_vol(i) * dz);
    end
    T_new = T_old + dT_dt * dt;
end

% 实现“追赶法”（Thomas Algorithm）
% 这是求解一维扩散/流动方程（产生的三对角矩阵）最高效的数值算法
function x = tridiag_solver(a, b, c, d)
    % 标准三对角求解器 (省略注释以节省空间)
    n = length(d);
    c_prime = zeros(n,1); d_prime = zeros(n,1); x = zeros(n,1);
    if b(1)==0, b(1)=1e-10; end 
    c_prime(1) = c(1)/b(1); d_prime(1) = d(1)/b(1);
    for i=2:n
        temp = b(i) - a(i)*c_prime(i-1);
        if temp==0, temp=1e-10; end
        c_prime(i) = c(i)/temp;
        d_prime(i) = (d(i) - a(i)*d_prime(i-1))/temp;
    end
    x(n) = d_prime(n);
    for i=n-1:-1:1
        x(i) = d_prime(i) - c_prime(i)*x(i+1);
    end
end