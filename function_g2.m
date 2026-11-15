function g2=function_g2(x,y)
% 速度分量 u2 的 Dirichlet 边界值（讲义第 49 页）
% 边界上直接取解析解 u2 = -(2/3) x y^3 + 2 - π sin(πx) 在该点的值。
% [修改7] 注释由原弹性方程齐次边界更新为 Stokes 解析解边界值；取值接口不变。

g2=exact_u_2D(x,y,2,0,0);

end
