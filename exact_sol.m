function sol_exact = exact_sol(x,y)
% Stokes 算例解析解向量 [u1; u2; p]（讲义第 50 页）
% [修改3] 原文件为弹性方程解析解，现替换为 Stokes 的 u1、u2、p。

u1=exact_u_2D(x,y,1,0,0);
u2=exact_u_2D(x,y,2,0,0);
p =exact_u_2D(x,y,3,0,0);
sol_exact=[u1;u2;p];
end
