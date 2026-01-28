function [NextDOM, NextSOM] = process_dom_sorption(DOM, SOM_Inert, Params, dt, Theta, SoilParams, layer_depth)
% PROCESS_DOM_SORPTION DOM 吸附/解吸过程 (C 和 N 耦合)
%
% 输入:
%   DOM: 结构体 {C, N} (向量，各层，单位 g/m2)
%   SOM_Inert: 结构体 {C, N} (向量，各层). 
%              注意: 在 Daisy 中吸附态通常被视为一种 "Sorb" 库或归入 SOM 库。
%              这里我们假设吸附到 SOM3 (Inert) 或单独的吸附库。
%   Params: Params.DOM (包含 K_ads, S_max, alpha, ModelType)
%   Theta: 当前时刻体积含水量 (向量, cm3/cm3)
%   SoilParams: 土壤参数结构体 (包含 BD, g/cm3)
%   layer_depth: 层厚 (cm)
% 输出:
%   NextDOM, NextSOM: 更新后的状态(g/m2)
% 2025-12-17 lc

    NextDOM = DOM;
    NextSOM = SOM_Inert;
    
    n_lay = length(DOM.C);
    
    % 获取参数
    alpha = Params.alpha; % 动力学速率 (1/h)
    
    for z = 1:n_lay
        sol_C_mass = DOM.C(z); % g/m2
        sol_N_mass = DOM.N(z);
        sorb_C_mass = SOM_Inert.C(z); % g/m2
        sorb_N_mass = SOM_Inert.N(z);

        % 0. 单位转换: 质量(g/m2) -> 浓度
        % 液相浓度 C_liq [g/cm3]
        % C_liq = (Mass_g_per_m2 * 1e-4) / (Depth_cm * Theta)
        th = max(1e-6, Theta(z));
        conc_sol_C = (sol_C_mass * 1e-4) / (layer_depth * th);
        
        % 1. 计算平衡吸附浓度 S_eq_conc [g solute / g soil]
        % 假设 S_max 单位是 g/g_soil, K_adsorption 单位是 cm3/g_solute (Langmuir standard)
       
        switch Params.SorptionModel
            case 'Langmuir'
                % S = (Smax * K * C) / (1 + K * C)
                % 注意: 这里输入的 conc_sol_C 必须是 g/cm3 (或 mg/L 需对应 K 单位)
                % 假设 K_adsorption 适配 g/cm3
                S_eq_conc = (Params.S_max * Params.K_adsorption * conc_sol_C) / ...
                       (1 + Params.K_adsorption * conc_sol_C);
                   
            case 'Linear'
                % S = Kd * C
                S_eq_conc = Params.K_dist * conc_sol_C;
                
            case 'Freundlich'
                % S = Kf * C^n
                S_eq_conc = Params.K_freundlich * (conc_sol_C ^ Params.n_freundlich);
                
            otherwise
                S_eq_conc = 0; % 默认无吸附
        end
        
       % 2. 将平衡浓度转回质量 [g/m2]
        % S_mass = S_conc * Mass_soil
        % Mass_soil (g/m2) = BD (g/cm3) * Depth (cm) * 1e4 (cm2/m2)
        mass_soil_per_m2 = SoilParams.BD(z) * layer_depth * 1e4;
        S_eq_mass = S_eq_conc * mass_soil_per_m2;
        
        % 3. 动力学传质: dS/dt = alpha * (S_eq - S)
        % 现在 S_eq_mass 和 sorb_C_mass 都是 g/m2，可以直接相减
        dC_sorp = alpha * (S_eq_mass - sorb_C_mass) * dt;
        
        % 4. 质量守恒限制
        if dC_sorp > 0 % 吸附 (Solution -> Solid)
            dC_sorp = min(dC_sorp, sol_C_mass);
            
            % 计算伴随的 N 吸附
            % 假设吸附的是 DOM 整体，携带其当前的 C/N 比
            if sol_C_mass > 1e-9
                current_CN = sol_C_mass / max(1e-9, sol_N_mass);
                dN_sorp = dC_sorp / current_CN;
            else
                dN_sorp = 0;
            end
            
        else % 解吸 (Solid -> Solution)
            dC_sorp = max(dC_sorp, -sorb_C_mass);
            
            % 计算伴随的 N 解吸
            % 假设解吸的是固相吸附物，携带其当前的 C/N 比
            if sorb_C_mass > 1e-9
                current_CN = sorb_C_mass / max(1e-9, sorb_N_mass);
                dN_sorp = dC_sorp / current_CN;
            else
                dN_sorp = 0;
            end
        end
        
        % 4. 更新状态
        NextDOM.C(z) = sol_C_mass - dC_sorp;
        NextDOM.N(z) = sol_N_mass - dN_sorp;
        
        NextSOM.C(z) = sorb_C_mass + dC_sorp;
        NextSOM.N(z) = sorb_N_mass + dN_sorp;
    end
end