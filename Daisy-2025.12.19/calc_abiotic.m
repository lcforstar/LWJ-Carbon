function Factors = calc_abiotic(T, h, clay, ClayModel)
% CALC_ABIOTIC 计算非生物环境因子 (温度, 水分, 粘土)
%
% 输入:
%   T: 温度 (摄氏度)
%   h: 水势 (cm 水柱, 负值)
%   clay: 粘土含量 (0-1, 例如 0.15 代表 15%)
%   ClayModel: 字符串, 选择粘土模型
%       - 'Standard': (默认) Daisy 经典模型 (ClayOMOld)
%       - 'Biomod': BIOMOD 项目模型 (ClayOMBiomod)
%       - 'None': 不考虑粘土效应
%
% 输出:
%   Factors.Turnover: 用于修正分解速率的因子
%   Factors.Maintenance: 用于修正微生物维持呼吸的因子
% T: 摄氏度, h: cm 水柱, clay: 0-1\
% 2025-12-16 lc

    
    % 1. 温度因子 f_T0 (Daisy standard)
    if T < 0
        f_T = 0;
    elseif T < 20
        f_T = 0.1 * T;
    elseif T < 37
        f_T = exp(0.47 - 0.027 * T + 0.00193 * T * T);
    elseif T < 60
        f_T_37 = exp(0.47 - 0.027 * 37.0 + 0.00193 * 37.0 * 37.0);
        f_T = f_T_37 * (1 - (T-37.0) / (60.0-37.0));
    else
        % T >= 60
        f_T = 0;
    end
    
    % 2. 水分因子 f_h (Daisy pF curve)
    if h >= 0
        f_W = 0.6; % 饱和抑制
    else
        pF = log10(-h);
        if pF <= 0
            f_W = 0.6;
        elseif pF <= 1.5
            f_W = 0.6 + (1.0 - 0.6) * pF / 1.5;
        elseif pF <= 2.5
            f_W = 1.0; % 最适
        elseif pF <= 6.5
            f_W = 1.0 - (pF - 2.5) / (6.5 - 2.5);
        else
            f_W = 0; % 极干
        end
    end
    
    Base = f_T * f_W;
    
    % 3. 粘土效应
    switch ClayModel
        case 'Standard'
            % 经典 Daisy: 粘土保护有机质(降低分解), 但增加微生物维持消耗
            % Turnover: 线性下降, 100% 粘土时降至 0.25
            f_clay_turnover = max(0, 1.0 - 0.75 * clay);
            % Maintenance: 线性上升, 粘土越多维持消耗越大
            f_clay_maint = 1.0 + clay;
            Factors.Turnover = Base * f_clay_turnover;
            Factors.Maintenance = Base * f_clay_maint;
            
        case 'Biomod'
            % BIOMOD 参数化: 使用 PLF (Piecewise Linear Function)
            % 逻辑: 粘土在 0-25% 之间强烈抑制分解，超过 25% 后效果饱和
            % 特点: 对维持呼吸也是抑制作用 (乘以相同因子)
            if clay <= 0.25
                % 从 (0.00, 1.0) 到 (0.25, 0.5) 线性下降
                % 斜率 = (0.5 - 1.0) / 0.25 = -2.0
                factor = 1.0 - 2.0 * clay;
            else
                % 大于 25% 粘土时保持 0.5
                factor = 0.5;
            end
            
            Factors.Turnover = Base * factor;
            Factors.Maintenance = Base * factor; % 注意这里也是降低
            
        case 'None'
            % 无粘土效应
            Factors.Turnover = Base;
            Factors.Maintenance = Base;
            
        otherwise
            error('未知的 ClayModel: %s', ClayModel);
    end
end