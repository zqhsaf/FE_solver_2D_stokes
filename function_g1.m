function g1=function_g1(x,y)
% 速度分量 u1 的 Dirichlet 边界值（讲义第 49 页）
% 边界上直接取解析解 u1 = x^2 y^2 + e^{-y} 在该点的值。
% [修改6] 注释由原弹性方程齐次边界更新为 Stokes 解析解边界值；取值接口不变。

g1=exact_u_2D(x,y,1,0,0);

end
