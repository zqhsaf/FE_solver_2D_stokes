function [solution,P,T,Pb_u,Tb_u,Pb_p,Tb_p]=FE_solver_2D_stokes(left,right,bottom,top,N1,N2)
% 二维稳态 Stokes 方程 Taylor-Hood 混合有限元求解器
%   速度 u：二次三角元（basis_type_u=202，6 个局部基）
%   压力 p：一次三角元（basis_type_p=201，3 个局部基）
% 弱形式（讲义第 11、15、20 页）：
%   a(u,v) + b(v,p) = (f,v),   b(u,q) = 0
%   a(u,v)=∫2νD(u):D(v) dxdy,  b(v,p)=-∫p ∇·v dxdy
% 算例区域 Ω=[0,1]×[-0.25,0]，解析解见 exact_u_2D.m。

basis_type_u=202;
basis_type_p=201;

%% 1. 生成网格与有限元节点
[P,T]=generate_PT(left,right,bottom,top,N1,N2);
[Pb_u,Tb_u]=generate_PbTb(left,right,bottom,top,N1,N2,basis_type_u);
[Pb_p,Tb_p]=generate_PbTb(left,right,bottom,top,N1,N2,basis_type_p);

number_of_local_basis_fun_u=6;
number_of_local_basis_fun_p=3;

%% 2. 生成边界节点（基于坐标，一次/二次元通用）
boundarynodes=generate_boundarynodes(Pb_u,left,right,bottom,top);

Nb =size(Pb_u,2);    % 速度有限元节点数（u1、u2 各 Nb 个自由度）
Nbp=size(Pb_p,2);    % 压力有限元节点数
% [修改9] 修复：原代码 matrix_size_u 之后多处误用未定义变量 matrix_size；
%         并按各耦合块的“检验基行数 × 试探基列数”补齐正确尺寸。
matrix_size_u =[Nb, Nb];    % 速度-速度块 A1~A4
matrix_size_up=[Nb, Nbp];   % 速度(检验)-压力(试探)块 A5、A6
matrix_size_pu=[Nbp, Nb];   % 压力(检验)-速度(试探)块 A7、A8
number_of_elements=size(T,2);

%% 3. 组装刚度矩阵（讲义第 29、42 页）
% A1=∫ν φj_x φi_x ; A2=∫ν φj_y φi_y
A1=assemble_matrix_2D('function_nu',matrix_size_u,number_of_elements,P,T, ...
    Tb_u,Tb_u,number_of_local_basis_fun_u,number_of_local_basis_fun_u, ...
    basis_type_u,1,0,basis_type_u,1,0);
A2=assemble_matrix_2D('function_nu',matrix_size_u,number_of_elements,P,T, ...
    Tb_u,Tb_u,number_of_local_basis_fun_u,number_of_local_basis_fun_u, ...
    basis_type_u,0,1,basis_type_u,0,1);
% A3=∫ν φj_x φi_y ; A4=∫ν φj_y φi_x
A3=assemble_matrix_2D('function_nu',matrix_size_u,number_of_elements,P,T, ...
    Tb_u,Tb_u,number_of_local_basis_fun_u,number_of_local_basis_fun_u, ...
    basis_type_u,1,0,basis_type_u,0,1);
A4=assemble_matrix_2D('function_nu',matrix_size_u,number_of_elements,P,T, ...
    Tb_u,Tb_u,number_of_local_basis_fun_u,number_of_local_basis_fun_u, ...
    basis_type_u,0,1,basis_type_u,1,0);
% A5=-∫ψj φi_x ; A6=-∫ψj φi_y （试探=压力，检验=速度）
A5=assemble_matrix_2D('function_negativeone',matrix_size_up,number_of_elements,P,T, ...
    Tb_p,Tb_u,number_of_local_basis_fun_p,number_of_local_basis_fun_u, ...
    basis_type_p,0,0,basis_type_u,1,0);
A6=assemble_matrix_2D('function_negativeone',matrix_size_up,number_of_elements,P,T, ...
    Tb_p,Tb_u,number_of_local_basis_fun_p,number_of_local_basis_fun_u, ...
    basis_type_p,0,0,basis_type_u,0,1);
% A7=-∫φj_x ψi ; A8=-∫φj_y ψi （试探=速度，检验=压力）
A7=assemble_matrix_2D('function_negativeone',matrix_size_pu,number_of_elements,P,T, ...
    Tb_u,Tb_p,number_of_local_basis_fun_u,number_of_local_basis_fun_p, ...
    basis_type_u,1,0,basis_type_p,0,0);
A8=assemble_matrix_2D('function_negativeone',matrix_size_pu,number_of_elements,P,T, ...
    Tb_u,Tb_p,number_of_local_basis_fun_u,number_of_local_basis_fun_p, ...
    basis_type_u,0,1,basis_type_p,0,0);

A_O=zeros(Nbp,Nbp);
% 整体分块矩阵（讲义第 29 页）：
% A = [2A1+A2, A3, A5; A4, 2A2+A1, A6; A7, A8, O]
A=[2*A1+A2, A3,        A5;
   A4,      2*A2+A1,   A6;
   A7,      A8,        A_O];

%% 4. 组装载荷向量 b = ∫ f·ψ（讲义第 31、44 页）
vector_size=[Nb,1];
b1=assemble_vector_2D('function_f1',vector_size,number_of_elements,P,T, ...
    Tb_u,number_of_local_basis_fun_u,basis_type_u,0,0);
b2=assemble_vector_2D('function_f2',vector_size,number_of_elements,P,T, ...
    Tb_u,number_of_local_basis_fun_u,basis_type_u,0,0);
b_O=zeros(Nbp,1);
b=[b1;b2;b_O];

%% 5. 处理速度 Dirichlet 边界条件（讲义第 36、45 页，只改前两块速度自由度）
% [修改10] 修复：原代码传入未定义的 Pb_test，应为速度节点矩阵 Pb_u。
[A,b]=treat_Dirichlet_boundary(A,b,boundarynodes,Pb_u,'function_g1','function_g2',-1);

%% 6. 压力唯一性处理：固定域内一个压力节点（讲义第 37、40 页方案 (1)）
% [修改11] 新增：压力只出现在散度项中，(u,p) 与 (u,p+c) 同为解，系统对压力
%         差一个常数零空间；需固定一个压力自由度消除奇异性。
%         解向量分段：u1(1:Nb)、u2(Nb+1:2Nb)、p(2Nb+1:2Nb+Nbp)。
tol=1e-10;
% 优先选取一个不在四条边界上的内部压力节点；若不存在则取第 1 个压力节点。
is_interior = abs(Pb_p(1,:)-left)>tol & abs(Pb_p(1,:)-right)>tol & ...
              abs(Pb_p(2,:)-bottom)>tol & abs(Pb_p(2,:)-top)>tol;
idx_interior=find(is_interior,1);
if isempty(idx_interior)
    k0=1;
else
    k0=idx_interior;
end
p_dof=2*Nb+k0;                       % 该压力节点对应的全局自由度编号
x0=Pb_p(1,k0); y0=Pb_p(2,k0);
A(p_dof,:)=0;
A(p_dof,p_dof)=1;
% 固定到解析压力在该点的值，使压力常数偏移最小，保证 L2/L∞ 按 O(h²) 收敛；
% （若置 0，压力会差一个不随 h 收敛的常数，L2/L∞ 无法达到二阶收敛。）
b(p_dof)=exact_u_2D(x0,y0,3,0,0);

%% 7. 求解整体线性方程组
solution=A\b;

end
