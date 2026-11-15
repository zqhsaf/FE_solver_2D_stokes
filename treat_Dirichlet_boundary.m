function [A,b] = treat_Dirichlet_boundary(A, b, boundarynodes, Pb, g1, g2, dirichlet_tag)
% 按“行置 1 法”处理向量值 Dirichlet 边界条件（Stokes 速度边界，讲义第 36、45 页）
% 每个速度网格节点有 2 个自由度：
%   第 k 个节点的 u1 自由度编号为 k，u2 自由度编号为 Nb+k（Nb=size(Pb,2)）
% 压力自由度不在本函数中处理（见 FE_solver_2D_stokes.m 第 6 步压力固定）。
% g1,g2：两个速度分量的边界值函数，可以是函数句柄或函数名字符串，
%        调用形式 g(x,y)，返回标量
% dirichlet_tag 为边界节点标记，默认 -1
% [修改12] 头部注释由弹性问题更新为 Stokes 问题；行置 1 算法逻辑保持不变。
%
% 历史诊断（早期弹性版本遗留，供参考）：原代码曾存在三处错误——
%   (1) 只调用一个 g(x,y)，而它返回的是 [g1;g2] 向量，赋给标量 b(i) 会报错；
%   (2) u2 自由度被赋成了 g1 的值（坐标也取错）；
%   (3) A(Nb+i,:)=0 写在循环里反复执行，且与 i 的处理混在一起。
%   现改为按节点 i 同时处理 k=i 与 Nb+i 两个自由度。

Nb=size(Pb,2);
if nargin < 7 || isempty(dirichlet_tag)
    dirichlet_tag = -1;
end

if ischar(g1) || isstring(g1)
    g1 = str2func(g1);
end
if ischar(g2) || isstring(g2)
    g2 = str2func(g2);
end

nbn = size(boundarynodes, 2);

for k = 1 : nbn

    if boundarynodes(1, k) == dirichlet_tag

        i = boundarynodes(2, k);        % 该边界节点的全局节点编号
        x = Pb(1, i);
        y = Pb(2, i);

        % ---- 第 1 个位移分量 u1：自由度 i ----
        A(i, :) = 0;                    % 边界自由度所在行清零
        A(i, i) = 1;                    % 对角元置 1
        b(i)    = g1(x, y);

        % ---- 第 2 个位移分量 u2：自由度 Nb+i ----
        j = Nb + i;
        A(j, :) = 0;
        A(j, j) = 1;
        b(j)    = g2(x, y);

    end
end

end
