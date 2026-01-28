function [C_new, J_flux] = process_transport(C_old, theta_old, theta_new, q_edge, dt, layer_depth, Params)
% PROCESS_TRANSPORT 溶质的对流与水动力弥散 (Convection-Dispersion)
% 基于 Daisy transport_Hansen.C 实现
%
% 输入:
%   C_old:      上一时刻溶质浓度 [g/cm3] (向量)
%   theta_old:  上一时刻体积含水量 [cm3/cm3]
%   theta_new:  当前时刻体积含水量
%   q_edge:     层间水通量 [cm/h], q_edge(i) 是第 i-1 层和第 i 层之间的流量。
%               长度应为 N+1 (包括上边界和下边界)
%               正值表示向下流动 (注意：Daisy 中有时定义向上为正，这里我们要统一方向，向下为正)
%   dt:         时间步长 [h]
%   layer_depth:每层厚度 [cm] (标量或向量)
%   Params:     参数结构体 (dispersivity (lambda)
%                           diffusion_coeff (D0)
%                           porosity (phi)
%
% 输出:
%   C_new:      更新后的浓度
%   J_flux:     各层间的溶质通量 [g/cm2/h]
% 2025-12-15 lc

    % 确定土层
    N = length(C_old); 
    %处理土层厚度。如果输入的厚度是一个数（比如10cm），就把它变成一个列表，每一层都是10cm；如果输入的就是列表，就直接用
    if isscalar(layer_depth), dz = ones(N,1) * layer_depth; else, dz = layer_depth; end
    
    % 确保土层界面 = 土层 +1
    if length(q_edge) ~= N+1
        error('q_edge 长度必须是 层数+1 (包含上下边界通量)');
    end

    % --- 1. 计算弥散系数 D [cm^2/h] ---
    % D = lambda * |v| + tau * D0
    % 对应 C++ 中的 calculation of D vector
    %创建空列表D，储存所有界面处的弥散系数，初始设置为0
    D = zeros(N+1, 1);
    lambda = Params.dispersivity; % 纵向弥散度 [cm]
    D0 = Params.diffusion_coeff;  % 分子扩散系数 [cm2/h]
    
    % 边界 D 处理: 
    % 上边界(Surface): 只有入渗，通常忽略地表之上的弥散，设为0或很小
    D(1) = 0; 

    %从第二个界面开始计算
    for i = 2:N+1
        % 界面处的特征含水量 (上下两层、新旧时刻值取平均，最后一层只取最后一层的值)
        if i <= N
            th_avg = (theta_old(i-1) + theta_old(i) + theta_new(i-1) + theta_new(i)) / 4.0;
        else
            th_avg = (theta_old(N) + theta_new(N)) / 2.0; % 下边界
        end
        th_avg = max(1e-6, th_avg); %防止含水量为0计算出错，这里设置最小值
        
        % 孔隙流速 v = q / theta
        v = q_edge(i) / th_avg;
        
        % 曲折度因子 (Millington-Quirk)，土越干，路越难走
        phi = Params.porosity;
        tau = (th_avg .^ (7/3)) ./ (phi^2);

        % 计算总弥散系数D
        % lambda * abs(v)水动力弥散 
        % tau * D0分子扩散
        D(i) = lambda * abs(v) + tau * D0; 
    end
    
    % --- 2. 构建三对角矩阵系数 ---
    % 求解方程: -a*C(i-1) + b*C(i) - c*C(i+1) = d
    % 使用全隐式 (Fully Implicit) 或 Crank-Nicolson
    % 这里复现 Hansen 的隐式逻辑

    % 准备四个空列表存储方程组系数
    a = zeros(N, 1);
    b = zeros(N, 1);
    c = zeros(N, 1);
    d = zeros(N, 1);
    
    % 上边界通量 (Flux input)
    J_top = q_edge(1) * 0; % 假设顶部入流浓度为 0 (除非有降雨带入 N)
    % 如果有降雨沉降氮，应在外部计算好 J_top
    if isfield(Params, 'C_rain') && q_edge(1) > 0
        J_top = q_edge(1) * Params.C_rain;
    end
    
    for i = 1:N
        % 几何参数
        dzi = dz(i); %获取当前土层厚度/深度
        % 计算当前层中心到上一层中心的距离
        if i > 1
            dz_up = (dz(i) + dz(i-1))/2; % 节点间距
        else
            dz_up = dz(i)/2; % 表层到节点
        end
         
        % 计算当前层中心到下一层中心的距离
        if i < N
            dz_dn = (dz(i) + dz(i+1))/2;
        else
            dz_dn = dz(i)/2;
        end
        
        % 流量权重 (Upwind scheme: alpha=1 if q>0, alpha=0 if q<0)
        % Daisy C++ 代码中: if (q < 0) alpha = 1.0 (因为它是 z 轴向上为正)
        % 我们这里 z 轴向下为正，所以 if q > 0 (下渗), alpha = 1.0 (取上游)
        alpha_up = (q_edge(i) > 0);
        alpha_dn = (q_edge(i+1) > 0);
        
        % 离散化项
        % 这一部分比较复杂，对应 C++ 中的 a[j], b[j], c[j], d[j] 赋值
        % 简化为标准有限差分格式：
        
        % Flux J = q*C - D*dC/dz
        % d(theta*C)/dt = -dJ/dz
        
        % 转换为线性方程组形式
        % V_i * (theta_new * C_new - theta_old * C_old) / dt = J_in - J_out
        
        % 这是一个非线性过程的线性化近似
        
        % 通量项系数 (Implicit, time = t_new)
        % J_up = q_up * (alpha*C_up + (1-alpha)*C_i) - D_up * (C_i - C_up)/dz_up
        % J_dn = q_dn * (alpha*C_i + (1-alpha)*C_dn) - D_dn * (C_dn - C_i)/dz_dn
        
        % 计算来自上一层的影响系数。如果是第 1 层，就用固定的 J_top
        if i == 1
            % 上边界处理 (Flux boundary)
            coef_Cm1 = 0; 
            flux_in_fixed = J_top; % 包含对流和弥散的固定通量输入
        else
        % 否则，既有水冲下来的（q*alpha），也有扩散下来的（D/dz）
            coef_Cm1 = q_edge(i) * alpha_up + D(i) / dz_up;
            flux_in_fixed = 0;
        end
        
        % 计算流出到上一层（或上一层反向流回来）的系数
        coef_C_in_up = q_edge(i) * (1 - alpha_up) - D(i) / dz_up;
        
        % 计算流出到下一层的系数
        coef_C_in_dn = q_edge(i+1) * alpha_dn + D(i+1) / dz_dn;
        
        % 计算来自下一层的影响系数。如果是最底层，假设水直接流走不回头（自由排水），所以系数修正一下
        if i == N
            % 下边界处理 (自由排水 dC/dz=0 -> J = q*C, 或者零浓度梯度)
            % 假设 q*C (纯对流流出，忽略底部弥散梯度)
            coef_Cp1 = 0;
            % 修正 coef_C_in_dn，因为没有 C(N+1)
            % J_out = q * C_N
            coef_C_in_dn = q_edge(i+1); 
        else
            coef_Cp1 = q_edge(i+1) * (1 - alpha_dn) - D(i+1) / dz_dn;
        end
        
        % 组装矩阵
        % LHS: theta_new * C_new / dt - (J_in_new - J_out_new) / dzi
        %      J_in_new = coef_Cm1 * C(i-1) + coef_C_in_up * C(i)
        %      J_out_new = coef_C_in_dn * C(i) + coef_Cp1 * C(i+1)
        
        % 方程: 
        % -coef_Cm1 * C(i-1) + (theta_new/dt*dzi - coef_C_in_up + coef_C_in_dn) * C(i) + coef_Cp1 * C(i+1) 
        % = theta_old * C_old / dt * dzi + flux_in_fixed
        
        % a 代表“上一层浓度”对方程的贡献
        % b 代表“本层浓度”的系数。包含水的变化（第一项）和进出通量的总和
        % c 代表“下一层浓度”对方程的贡献
        % d 是方程右边的已知数。主要是“原来有多少货”（旧浓度）加上“外部额外加的货”（flux_in_fixed）

        if i > 1
            a(i) = -coef_Cm1; 
        end
        
        b(i) = (theta_new(i) * dzi / dt) - coef_C_in_up + coef_C_in_dn; 
        
        if i < N
            c(i) = coef_Cp1; % 注意符号，通常移到左边是正的，因为 J_out 是减去 
        end
        
        d(i) = (theta_old(i) * C_old(i) * dzi / dt) + flux_in_fixed;
    end
    
    % --- 3. 调用“三对角矩阵求解器”求解 ---
    C_new = tridiag(a, b, c, d);
    
    % --- 4. 计算通量 (用于输出) ---
    % 初始化通量列表。第 1 个通量就是最上面进来的量
    J_flux = zeros(N+1, 1);
    J_flux(1) = J_top;
    % 算浓度梯度（浓度差 ÷ 距离）。这是为了算扩散
    for i = 2:N+1
        if i <= N
            c_grad = (C_new(i) - C_new(i-1)) / ((dz(i)+dz(i-1))/2);
            % 算随水流动的浓度
            % 如果水向下流（>0），这一流水的浓度等于上一层的浓度（C_new(i-1)
            % 如果水向上流，浓度等于本层的浓度
            c_adv = (q_edge(i)>0) * C_new(i-1) + (q_edge(i)<=0) * C_new(i); 
        else
            c_grad = 0; % 底部边界  
            c_adv = C_new(N);
        end
        % 计算总通量
        J_flux(i) = q_edge(i) * c_adv - D(i) * c_grad;
    end
end