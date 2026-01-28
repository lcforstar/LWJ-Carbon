function [NextState, Fluxes] = daisy_step(CurrentState, Params, Env, dt)
% 核心步长函数：包含 运移(Transport)、吸附(Sorption)、生物周转(Turnover)
% 单位说明: 
%   - 内部状态 (State) 存储为 g/m2 (质量)
%   - 运移过程 (Transport) 需要 g/cm3 (浓度)
%   - 代码中包含 mass2conc 和 conc2mass 转换

% 2025-12-18 lc

    NextState = CurrentState;
    n_layers = size(CurrentState.SOM.C, 1);
    
    % 获取库的数量
    n_smb = size(CurrentState.SMB.C, 2);
    n_som = size(CurrentState.SOM.C, 2);
    n_aom = size(CurrentState.AOM.C, 2);

    layer_depth = 10; % cm (假设均匀分层，需与 main 保持一致)
    dz_vec = ones(n_layers, 1) * layer_depth; % 向量化深度

    Fluxes.CO2_Total = zeros(n_layers, 1);
    Fluxes.Leaching.NO3 = 0;
    Fluxes.Leaching.DOM = 0;

    % A. 溶质运移 (Solute Transport) 
    % 使用对流-弥散方程 (CDE) 模拟物质随水流在土层间的移动
    % 注意: 只有溶解态物质移动 (DOM, NO3, NH4)
    % DAISY C++ 内核使用 CGS (cm, g, h)。
    % 转换: 1 g/m2 = 1e-4 g/cm2.
    
    % 辅助函数定义
    % M: g/m2, theta: cm3/cm3, dz: cm -> C: g/cm3
    mass2conc = @(M, th, dz) (M .* 1e-4) ./ (dz .* max(1e-6, th));
    % C: g/cm3, theta: cm3/cm3, dz: cm -> M: g/m2
    conc2mass = @(C, th, dz) C .* (dz .* max(1e-6, th)) .* 1e4;
   % 1. DOM 运移 (C 和 N 分别运移)
    % (1) 转为浓度
    Conc_DOM_C = mass2conc(CurrentState.DOM.C, Env.Theta_old, dz_vec);
    Conc_DOM_N = mass2conc(CurrentState.DOM.N, Env.Theta_old, dz_vec);
    
    % (2) 执行运移
    [Conc_DOM_C_new, J_DOM_C] = process_transport(...
        Conc_DOM_C, Env.Theta_old, Env.Theta_new, Env.q_edge, ...
        dt, layer_depth, Params.Transport);
        
    [Conc_DOM_N_new, J_DOM_N] = process_transport(...
        Conc_DOM_N, Env.Theta_old, Env.Theta_new, Env.q_edge, ...
        dt, layer_depth, Params.Transport);
    
    % (3) 转回质量
    NextState.DOM.C = conc2mass(Conc_DOM_C_new, Env.Theta_new, dz_vec);
    NextState.DOM.N = conc2mass(Conc_DOM_N_new, Env.Theta_new, dz_vec);
        
    % 2. 矿质氮运移 (NO3 移动性强，NH4 会被吸附移动性弱)
    % (1) 转为浓度
    Conc_NO3 = mass2conc(CurrentState.Soil.NO3, Env.Theta_old, dz_vec);
    
    % (2) 执行运移
    [Conc_NO3_new, J_NO3] = process_transport(...
        Conc_NO3, Env.Theta_old, Env.Theta_new, Env.q_edge, ...
        dt, layer_depth, Params.Transport);
        
    % (3) 转回质量
    NextState.Soil.NO3 = conc2mass(Conc_NO3_new, Env.Theta_new, dz_vec);
        
    % NH4 包含线性吸附 (Retardation Factor R = 1 + rho*Kd/theta)
    % 简化起见，这里假设 NH4 也随水移动，或者你可以设置 Retardation
    % (1) 转换为浓度 (g/m2 -> g/cm3)
    Conc_NH4 = mass2conc(CurrentState.Soil.NH4, Env.Theta_old, dz_vec);
    
    % (2) 执行运移
    [Conc_NH4_new, J_NH4] = process_transport(...
        Conc_NH4, Env.Theta_old, Env.Theta_new, Env.q_edge, ...
        dt, layer_depth, Params.Transport);
        
    % (3) 转回质量 (g/cm3 -> g/m2)
    NextState.Soil.NH4 = conc2mass(Conc_NH4_new, Env.Theta_new, dz_vec);

    % 记录淋溶流失量 (取最底部的通量)
    % process_transport 返回的 J 是 g/cm2/h，需转为 g/m2/h
    Fluxes.Leaching.NO3 = J_NO3(end) * 1e4;
    Fluxes.Leaching.NH4 = J_NH4(end) * 1e4;
    Fluxes.Leaching.DOM_C = J_DOM_C(end) * 1e4;
    Fluxes.Leaching.DOM_N = J_DOM_N(end) * 1e4;
    
    % B. 物理/化学过程 (Physics & Chemistry)
    
    % 1. AOM 生物混合 (Bioturbation) - 仅针对固体有机质
    [NextState.AOM.C, NextState.AOM.N] = process_bioinc(NextState.AOM.C, NextState.AOM.N, Env, dt);
    
    % 2. DOM 吸附/解吸 (Sorption) - C 与 N 耦合
    % 构造临时结构体以匹配函数接口
    DomPool.C = NextState.DOM.C; 
    DomPool.N = NextState.DOM.N;
    
    % 假设吸附态存储在 DOM.C_Sorb / N_Sorb 中 (也可映射到 SOM3)
    SorbPool.C = NextState.DOM.C_Sorb; 
    SorbPool.N = NextState.DOM.N_Sorb;
    
    % 调用吸附函数 
    [NewDom, NewSorb] = process_dom_sorption(DomPool, SorbPool, Params.DOM, dt, ...
                                             Env.Theta_new, Params.Soil, layer_depth);
    
    % 更新状态
    NextState.DOM.C = NewDom.C; NextState.DOM.N = NewDom.N;
    NextState.DOM.C_Sorb = NewSorb.C; NextState.DOM.N_Sorb = NewSorb.N;

    % C. 生化反应 (Biological Turnover)
    % 逐层计算微生物分解、矿化/固定
    for z = 1:n_layers
        % 1. 计算环境因子
        Abiotic = calc_abiotic(Env.Temp(z), Env.h(z), Env.Clay(z), Params.ClayModel);
        
 
        %【第一步】计算潜在分解 (Potential Decay)
        % AOM, SOM, SMB 的潜在分解
        Pot_AOM = get_potential(CurrentState.AOM.C(z,:), CurrentState.AOM.N(z,:), Params.AOM, Abiotic, dt);
        Pot_SOM = get_potential(CurrentState.SOM.C(z,:), CurrentState.SOM.N(z,:), Params.SOM, Abiotic, dt);
        Pot_SMB = get_potential_smb(CurrentState.SMB.C(z,:), CurrentState.SMB.N(z,:), Params.SMB, Abiotic, dt);
        Pot_DOM = get_potential_dom(NextState.DOM.C(z), NextState.DOM.N(z), Params.DOM, Abiotic, dt);
        
        %【第二步】氮平衡检查与减速因子 (Retardation)
        % 总 N 供应 (Supply) = 所有有机库分解释放的 N + 土壤矿质 N
        N_Supply_Organic = sum(Pot_AOM.N_out) + sum(Pot_SOM.N_out) + ...
                           sum(Pot_SMB.N_out) + sum(Pot_DOM.N_out);                 
        N_Available_Mineral = CurrentState.Soil.NH4(z) + CurrentState.Soil.NO3(z);
        
        % 总 N 需求 (Demand) = 微生物同化碳生成新生物量所需的 N
        % 注意: DOM 分解也会产生微生物生长需求
        N_Demand = calc_N_demand(Pot_AOM.C_out, Params.AOM, Params) + ...
                   calc_N_demand(Pot_SOM.C_out, Params.SOM, Params) + ...
                   calc_N_demand(Pot_SMB.C_out_turnover, Params.SMB, Params) + ...
                   calc_N_demand(Pot_DOM.C_out, Params.DOM, Params);
        
        % 计算减速因子 f_nit (0 ~ 1)
        f_nit = 1.0;
        if N_Demand > (N_Supply_Organic + N_Available_Mineral) + 1e-10
            f_nit = (N_Supply_Organic + N_Available_Mineral) / N_Demand;
            f_nit = max(0.01, f_nit); % 设置下限，防止完全停止
        end
        
        %【第三步】执行实际分解 (应用 f_nit)
        % 初始化本层增量缓冲
        dC_SMB = zeros(1, n_smb); dN_SMB = zeros(1, n_smb);
        dC_SOM = zeros(1, n_som); dN_SOM = zeros(1, n_som);
        dC_AOM = zeros(1, n_aom); dN_AOM = zeros(1, n_aom);
        dC_DOM = 0; dN_DOM = 0;   % DOM 是标量(每层)
        
        N_Min_Total = 0;
        Flux_CO2_z = 0;
        
        % 执行 AOM
        [dS, dT, co2, nmin] = process_aom(CurrentState.AOM.C(z,:), CurrentState.AOM.N(z,:), ...
                                          Params.AOM, Params, Abiotic, dt, f_nit);
        dC_AOM = dC_AOM + dS.dC; dN_AOM = dN_AOM + dS.dN;
        dC_SMB = dC_SMB + dT.SMB_dC; dN_SMB = dN_SMB + dT.SMB_dN;
        dC_SOM = dC_SOM + dT.SOM_dC; dN_SOM = dN_SOM + dT.SOM_dN;
        Flux_CO2_z = Flux_CO2_z + co2;
        N_Min_Total = N_Min_Total + nmin;
        
        % 执行 SOM
        [dS, dT, co2, nmin] = process_som(CurrentState.SOM.C(z,:), CurrentState.SOM.N(z,:), ...
                                          Params.SOM, Params, Abiotic, dt, f_nit);
        dC_SOM = dC_SOM + dS.dC; dN_SOM = dN_SOM + dS.dN;
        dC_SMB = dC_SMB + dT.SMB_dC; dN_SMB = dN_SMB + dT.SMB_dN;
        dC_SOM = dC_SOM + dT.SOM_dC; dN_SOM = dN_SOM + dT.SOM_dN;
        Flux_CO2_z = Flux_CO2_z + co2;
        N_Min_Total = N_Min_Total + nmin;
        
        % 执行 SMB
        [dS, dT, co2, nmin] = process_smb(CurrentState.SMB.C(z,:), CurrentState.SMB.N(z,:), ...
                                          Params.SMB, Params, Abiotic, dt, f_nit);
        dC_SMB = dC_SMB + dS.dC; dN_SMB = dN_SMB + dS.dN;
        dC_SMB = dC_SMB + dT.SMB_dC; dN_SMB = dN_SMB + dT.SMB_dN;
        dC_SOM = dC_SOM + dT.SOM_dC; dN_SOM = dN_SOM + dT.SOM_dN;
        Flux_CO2_z = Flux_CO2_z + co2;
        N_Min_Total = N_Min_Total + nmin;
        
       % 执行 DOM
        [dD, dS_D, co2, nmin] = process_dom_turnover(NextState.DOM.C(z), NextState.DOM.N(z), ...
                                                     Params.DOM, Params, Abiotic, dt, f_nit);
        dC_DOM = dC_DOM + dD.C; dN_DOM = dN_DOM + dD.N;
        dC_SMB = dC_SMB + dS_D.C; dN_SMB = dN_SMB + dS_D.N;
        Flux_CO2_z = Flux_CO2_z + co2;
        N_Min_Total = N_Min_Total + nmin;

        %【第四步】更新矿质氮 (Mineral N Update)
        if N_Min_Total > 0
            % 情况 1：净矿化
            NextState.Soil.NH4(z) = NextState.Soil.NH4(z) + N_Min_Total;
        else
            % 情况 2：净固定
            demand = -N_Min_Total;
            % 优先消耗 NH4
            if NextState.Soil.NH4(z) >= demand
                NextState.Soil.NH4(z) = NextState.Soil.NH4(z) - demand;
            else
                rem = demand - NextState.Soil.NH4(z);
                NextState.Soil.NH4(z) = 0;
                % 其次消耗 NO3
                if NextState.Soil.NO3(z) >= rem
                    NextState.Soil.NO3(z) = NextState.Soil.NO3(z) - rem;
                else
                    NextState.Soil.NO3(z) = 0;
                end
            end
        end
        
        % F. 应用所有库的增量 (Apply Deltas)
        Fluxes.CO2_Total(z) = Flux_CO2_z;
        
        NextState.AOM.C(z,:) = max(0, NextState.AOM.C(z,:) + dC_AOM);
        NextState.AOM.N(z,:) = max(0, NextState.AOM.N(z,:) + dN_AOM);
        
        NextState.SMB.C(z,:) = max(0, NextState.SMB.C(z,:) + dC_SMB);
        NextState.SMB.N(z,:) = max(0, NextState.SMB.N(z,:) + dN_SMB);
        
        NextState.SOM.C(z,:) = max(0, NextState.SOM.C(z,:) + dC_SOM);
        NextState.SOM.N(z,:) = max(0, NextState.SOM.N(z,:) + dN_SOM);
        
        NextState.DOM.C(z)   = max(0, NextState.DOM.C(z) + dC_DOM);
        NextState.DOM.N(z)   = max(0, NextState.DOM.N(z) + dN_DOM);
    end
end

% 辅助函数

function Pot = get_potential(C_pool, N_pool, PoolParams, Abiotic, dt)
    rate = min(PoolParams.k_turnover .* Abiotic.Turnover * dt, 1.0);
    Pot.C_out = C_pool .* rate;
    Pot.N_out = Pot.C_out ./ (C_pool ./ max(1e-9, N_pool)); 
    Pot.N_out(C_pool <= 1e-9) = 0;
end

function Pot = get_potential_smb(C_pool, N_pool, PoolParams, Abiotic, dt)
   rate_m = min(PoolParams.k_maint .* Abiotic.Maintenance * dt, 1.0);
    rate_t = min(PoolParams.k_turnover .* Abiotic.Turnover * dt, 1.0);
    C_maint = C_pool .* rate_m;
    Pot.C_out_turnover = max(0, C_pool - C_maint) .* rate_t;
    CN = C_pool ./ max(1e-9, N_pool);
    Pot.N_out = (C_maint + Pot.C_out_turnover) ./ CN; 
    Pot.N_out(C_pool<=1e-9)=0;
end

function Pot = get_potential_dom(C_dom, N_dom, PomParams, Abiotic, dt)
   rate = min(PomParams.k_turnover * Abiotic.Turnover * dt, 1.0);
    Pot.C_out = C_dom * rate;
    Pot.N_out = Pot.C_out / (C_dom / max(1e-9, N_dom)); 
    if C_dom<=1e-9, Pot.N_out=0; end
end

function N_req = calc_N_demand(C_out_vec, SourceParams, GlobalParams)
% CALC_N_DEMAND 计算碳分解后，生成微生物生物量所需的氮量
    N_req = 0;
    
    % 1. 提取参数
    Eff_all = SourceParams.efficiency;
    Fracs_all = SourceParams.fractions;
    
    % 2. 遍历每一个分解释放源 (例如 AOM1, AOM2...)
    for i = 1:length(C_out_vec)
        C_out = C_out_vec(i);
        
        % 如果该库没有分解，直接跳过
        if C_out <= 0, continue; end
        
        % A. 获取当前库的利用效率 (e)
        if isscalar(Eff_all)
            e = Eff_all; % 如果参数只有一个值，应用于所有库
        else
            e = Eff_all(i); % 如果参数是向量，取对应库的值
        end
        
        % B. 获取当前库的分配向量 (f_row)
        % 判断 Fracs_all 是矩阵(多行) 还是 向量(单行)
        if size(Fracs_all, 1) > 1
            % 如果是多行矩阵 (如 AOM, SMB)，取第 i 行
            f_row = Fracs_all(i, :);
        else
            % 如果是单行向量 (如 DOM)，直接使用
            f_row = Fracs_all;
        end
        
        % 3. 计算同化进入生物量的碳 (Assimilated Carbon)
        C_assim = C_out * e; % 同化碳量 = 分解量 * 利用效率
        
        % 4. 遍历所有去向 (SMB1, SMB2, SOM1...) 计算氮需求
        for t = 1:length(f_row)
            fraction = f_row(t);
            
            if fraction > 0
                % 根据 GlobalParams 确定目标库的 C/N 比
                % 索引对应关系: 1-2->SMB, 3-5->SOM
                if t <= 2 
                    target_CN = GlobalParams.SMB.CN(t);
                elseif t <= 5 
                    idx_som = t - 2;
                    target_CN = GlobalParams.SOM.CN(idx_som);
                else
                    target_CN = 10; % 默认保护值
                end
                
                % 累加氮需求: (进入该库的碳) / (该库的C/N)
                N_req = N_req + (C_assim * fraction) / target_CN; % N需求 = 同化碳 / 目标库C/N
            end
        end
    end
end