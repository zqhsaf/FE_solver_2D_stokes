function f1=function_f1(x,y)
% Stokes 问题右端项 f 的第一个分量 f1（讲义第 50 页）
%   f1 = -2ν x^2 - 2ν y^2 - ν e^{-y} + π^2 cos(πx) cos(2πy)
% [修改4] 原文件是弹性方程右端项，且引用了未定义的 lam/mu/function_lamda/
%         function_mu，无法运行；现按讲义 Stokes 算例整体替换。

nu = function_nu(x,y);

f1 = -2*nu*x.^2 - 2*nu*y.^2 - nu*exp(-y) ...
     + pi^2*cos(pi*x).*cos(2*pi*y);

end
