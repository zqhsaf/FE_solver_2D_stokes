function f2=function_f2(x,y)
% Stokes 问题右端项 f 的第二个分量 f2（讲义第 50 页）
%   f2 = 4ν x y - ν π^3 sin(πx) + 2π(2 - π sin(πx)) sin(2πy)
% [修改5] 原文件是弹性方程右端项，引用未定义的 lam/mu；现按讲义 Stokes 算例替换。

nu = function_nu(x,y);

f2 = 4*nu*x.*y - nu*pi^3*sin(pi*x) ...
     + 2*pi*(2 - pi*sin(pi*x)).*sin(2*pi*y);

end
