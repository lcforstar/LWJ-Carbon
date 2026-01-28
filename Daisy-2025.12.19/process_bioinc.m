function [C_AOM, N_AOM] = process_bioinc(C_AOM, N_AOM, Env, dt)
% 对应 bioincorporation.C
% 将第一层的 AOM 移动到更深层
% 2025-12-15 lc
    
    n_lay = size(C_AOM, 1);
    if n_lay < 2, return; end
    
    % 环境阈值
    if Env.Temp(1) < 4 || Env.h(1) < -1000, return; end
    
    % 混合速率 (Q10)
    k_bio = 0.0002; % 基础速率
    f_temp = 2^((Env.Temp(1)-10)/10);
    rate = k_bio * f_temp * dt;
    rate = min(rate, 0.5);
    
    % 搬运
    flux_C = C_AOM(1, :) * rate;
    flux_N = N_AOM(1, :) * rate;
    
    C_AOM(1, :) = C_AOM(1, :) - flux_C;
    N_AOM(1, :) = N_AOM(1, :) - flux_N;
    
    % 分配到 2-3 层
    target_z = 2:min(3, n_lay);
    n_t = length(target_z);
    for z = target_z
        C_AOM(z, :) = C_AOM(z, :) + flux_C / n_t;
        N_AOM(z, :) = N_AOM(z, :) + flux_N / n_t;
    end
end