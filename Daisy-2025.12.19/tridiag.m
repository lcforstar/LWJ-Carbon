function x = tridiag(a, b, c, d)
% TRIDIAG 求解三对角线性方程组 (Thomas Algorithm)
% 对应 Daisy C++ 中的 tridia 函数
% 方程形式: a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i)
%
% 输入:
%   a: 下对角线 (长度 N, a(1) 不使用)
%   b: 主对角线 (长度 N)
%   c: 上对角线 (长度 N, c(N) 不使用)
%   d: 常数项 (长度 N)
%
% 输出:
%   x: 解向量
% 2025-12-15 lc

    n = length(d);
    x = zeros(n, 1);
    
    % 消元过程 (Forward elimination)
    % 修改系数，使 a 变为 0
    for i = 2:n
        factor = a(i) / b(i-1);
        b(i) = b(i) - factor * c(i-1);
        d(i) = d(i) - factor * d(i-1);
    end
    
    % 回代过程 (Back substitution)
    x(n) = d(n) / b(n);
    for i = n-1:-1:1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i);
    end
end