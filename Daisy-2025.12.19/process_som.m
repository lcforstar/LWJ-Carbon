function [dS, dT, CO2, N_min] = process_som(C, N, Params, Global, Abiotic, dt, f_nit)
% 输入：
    % C,N: 当前 SOM 库现有碳氮质量
    % Params: SOM 的参数（分解速率 k，利用效率 e，分配比例）
    % Global：全局参数，获取目标库的碳氮比
    % Abiotic：环境因子（温度、水分、粘土含量）
    % dt：时间步长（一般一小时）
    % f_nit：氮限制因子（0~1）
% 输出：
    % dS：源库（SOM）的减少量
    % dT：目标库（SMB、SOM）的增加量
    % CO2:呼吸释放的CO2
    % N_min：净矿化量

% 2025-12-15 lc
    n_smb = 2; n_som = 3;
    dS.dC = zeros(1, n_som); dS.dN = zeros(1, n_som);
    dT.SMB_dC = zeros(1, n_smb); dT.SMB_dN = zeros(1, n_smb);
    dT.SOM_dC = zeros(1, n_som); dT.SOM_dN = zeros(1, n_som);
    
    CO2 = 0; N_min = 0;
    
    for i = 1:n_som
        if C(i) <= 1e-9, continue; end
        
        CN = C(i) / max(1e-9, N(i));
        
        rate = Params.k_turnover(i) * Abiotic.Turnover * f_nit;
        rate = min(rate * dt, 1.0);
        
        C_decay = C(i) * rate;
        N_decay = C_decay / CN;
        
        dS.dC(i) = -C_decay;
        dS.dN(i) = -N_decay;
        
        eff = Params.efficiency(i);
        C_assim = C_decay * eff;
        CO2 = CO2 + C_decay * (1 - eff);
        
        N_req_tot = 0;
        fracs = Params.fractions(i, :);
        
        for t = 1:5
            if fracs(t) > 0
                C_in = C_assim * fracs(t);
                if t <= 2 % SMB
                    target_CN = Global.SMB.CN(t);
                    N_req = C_in / target_CN;
                    dT.SMB_dC(t) = dT.SMB_dC(t) + C_in;
                    dT.SMB_dN(t) = dT.SMB_dN(t) + N_req;
                else % SOM (Index 3->1, 4->2, 5->3)
                    idx = t - 2;
                    target_CN = Global.SOM.CN(idx);
                    N_req = C_in / target_CN;
                    dT.SOM_dC(idx) = dT.SOM_dC(idx) + C_in;
                    dT.SOM_dN(idx) = dT.SOM_dN(idx) + N_req;
                end
                N_req_tot = N_req_tot + N_req;
            end
        end
        N_min = N_min + (N_decay - N_req_tot);
    end
end