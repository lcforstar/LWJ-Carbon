function [dDOM, dSMB, CO2, N_min] = process_dom_turnover(DC, DN, P, Params, Abio, dt, f_nit)
% PROCESS_DOM_TURNOVER DOM 的生物分解过程
% 输入:
    %   DC, DN: DOM 的碳、氮含量
    %   P: Params.DOM (分解速率 k，利用效率 e，分配比例)
    %   Params: 全局参数 (用于获取 SMB C/N)
    %   Abio: 环境因子
    %   dt: 时间步长
    %   f_nit: 氮限制因子 (0-1) 
% 输出：
    % dS：源库（DOM）的减少量
    % dT：目标库（SMB）的增加量
    % CO2:呼吸释放的CO2
    % N_min：净矿化量
% 2025-12-15 lc

    % 1. 计算分解速率 (受 f_nit 限制)
    rate = P.k_turnover * Abio.Turnover * f_nit; 
    rate = min(rate * dt, 1.0);
    
    % 2. 计算流失量
    C_decay = DC * rate;
    % 保持当前 C/N 比流失
    if DC > 1e-9
        N_decay = C_decay / (DC / max(1e-9, DN));
    else
        N_decay = 0;
    end
    
    % DOM 库的净改变 (流出)
    dDOM.C = -C_decay; 
    dDOM.N = -N_decay;
    
    % 3. 碳的分配 (呼吸 + 同化)
    CO2 = C_decay * (1 - P.efficiency);
    C_assim = C_decay * P.efficiency;
    
    % 4. 分配到 SMB (微生物生物量)
    % 处理 fractions 可能是标量或向量的情况
    n_smb = length(Params.SMB.CN);
    dSMB.C = zeros(1, n_smb); 
    dSMB.N = zeros(1, n_smb);
    
    N_req_tot = 0;
    
    % 确保分配矩阵为一维向量（因为DOM 不分层级，一对一分配）
    fracs = P.fractions;
    if size(fracs, 1) > 1, fracs = fracs(1, :); end % 防御性代码
    
    % 通常 DOM 主要流向SMB2（Fast）
    for i = 1:min(length(fracs), n_smb)
        if fracs(i) > 0
            C_in = C_assim * fracs(i);
            
            % 目标 SMB 库的 N 需求
            target_CN = Params.SMB.CN(i);
            N_req = C_in / target_CN;
            
            dSMB.C(i) = C_in; 
            dSMB.N(i) = N_req;
            
            N_req_tot = N_req_tot + N_req;
        end
    end
    
    % 5. 计算矿化/固定 (释放出的N - 微生物生长需要的N)
    N_min = N_decay - N_req_tot;
end