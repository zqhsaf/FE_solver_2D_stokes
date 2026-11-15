function nu = function_nu(x,y)
% Stokes 问题的动力粘度 ν（常数）
% [修改1] 新增：原求解器 assemble_matrix_2D 中以系数函数名 'function_nu' 调用，
%         但目录下缺少该文件；讲义第 4、42 页中 a(u,v)=∫2νD(u):D(v)，
%         本算例取常数粘度 ν=1（与讲义数值例 f1、f2 的系数一致）。
% 输入 x,y 仅为保持系数函数接口统一（逐元素可用），此处不依赖坐标。

nu = 1;

end
