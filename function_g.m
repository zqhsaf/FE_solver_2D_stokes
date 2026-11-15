function g=function_g(x,y)
% Stokes 问题 Dirichlet 边界上的速度向量 g=[g1;g2]（供检查/画图使用）。
% 求解时由 function_g1、function_g2 分别给定两个分量。
% [修改8] 注释更新为 Stokes 算例；内部仍分别调用 g1、g2。

g=[function_g1(x,y); function_g2(x,y)];

end
