function [dS, dT, CO2, N_min] = process_aom(C, N, Params, Global, Abiotic, dt, f_nit)
% AOM 逻辑与 SOM 几乎一致，只是来源不同
% 输入：
    % C,N: 当前 AOM 库现有碳氮质量
    % Params: AOM 的参数（分解速率 k，利用效率 e，分配比例）
    % Global：全局参数，获取目标库的碳氮比
    % Abiotic：环境因子（温度、水分、粘土含量）
    % dt：时间步长（一般一小时）
    % f_nit：氮限制因子（0~1）
% 输出：
    % dS：源库（AOM）的减少量
    % dT：目标库（SMB、SOM）的增加量
    % CO2:呼吸释放的CO2
    % N_min：净矿化量

% 2025-12-18 lc

    n_aom = 2; n_smb = 2; n_som = 3;
    % 初始化源库和目标库变化量结构体，创捷全零行向量，记录数据
    dS.dC = zeros(1, n_aom); dS.dN = zeros(1, n_aom);
    dT.SMB_dC = zeros(1, n_smb); dT.SMB_dN = zeros(1, n_smb);
    dT.SOM_dC = zeros(1, n_som); dT.SOM_dN = zeros(1, n_som);
    
    %记录本步长产生的总呼吸和总净矿化量
    CO2 = 0; N_min = 0;
    
    for i = 1:n_aom
        % 零值保护，库里没东西或极小值则跳过
        if C(i) <= 1e-9, continue; end
        CN = C(i) / max(1e-9, N(i)); % 计算当前库碳氮比
        
        rate = Params.k_turnover(i) * Abiotic.Turnover * f_nit; % 计算实际分解速率
        rate = min(rate * dt, 1.0); % 时间步长约束，防止分解量大于库总量
        
        C_decay = C(i) * rate; % 时间步长内分解的 C
        N_decay = C_decay / CN; % 伴随分解释放的 N
        
        % 记录 AOM 库的亏随
        dS.dC(i) = -C_decay;
        dS.dN(i) = -N_decay;
        
        eff = Params.efficiency(i); % 微生物利用效率 i = 0.4 表示40% 同化，60% 呼吸
        C_assim = C_decay * eff; % 计算同化碳量
        CO2 = CO2 + C_decay * (1 - eff); % 计算并累加呼吸碳量
        
        N_req_tot = 0; % 初始化 AOM 库分解过程的总氮需求
        fracs = Params.fractions(i, :); % 读取分配矩阵
        
        for t = 1:5
            if fracs(t) > 0
                C_in = C_assim * fracs(t);
                if t <= 2
                    target_CN = Global.SMB.CN(t);
                    dT.SMB_dC(t) = dT.SMB_dC(t) + C_in;
                    dT.SMB_dN(t) = dT.SMB_dN(t) + C_in/target_CN;
                    N_req_tot = N_req_tot + C_in/target_CN;
                else
                    idx = t - 2;
                    target_CN = Global.SOM.CN(idx);
                    dT.SOM_dC(idx) = dT.SOM_dC(idx) + C_in;
                    dT.SOM_dN(idx) = dT.SOM_dN(idx) + C_in/target_CN;
                    N_req_tot = N_req_tot + C_in/target_CN;
                end
            end
        end
        N_min = N_min + (N_decay - N_req_tot);
    end
end