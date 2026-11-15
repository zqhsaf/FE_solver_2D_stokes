function [u_Linf,u_L2,u_H1,p_Linf,p_L2,p_H1]=compute_error_2D(solution,P,T,Pb_u,Tb_u,basis_u,Pb_p,Tb_p,basis_p)
% 计算 Stokes 有限元数值解相对解析解的误差（讲义第 46~48 页定义）
%
% 输入：
%   solution ：(2*Nb+Nbp) 维列向量，分段为 u1(1:Nb)、u2(Nb+1:2Nb)、p(2Nb+1:end)
%   Pb_u/Tb_u：速度有限元节点坐标与连接表（二次元 202）
%   Pb_p/Tb_p：压力有限元节点坐标与连接表（一次元 201）
%
% 输出（速度 u 对应讲义第 52 页表，压力 p 对应第 53 页表）：
%   u_Linf = max(||u1-u1h||_inf, ||u2-u2h||_inf)
%   u_L2   = sqrt(||u1-u1h||_0^2 + ||u2-u2h||_0^2)
%   u_H1   = sqrt(|u1-u1h|_1^2  + |u2-u2h|_1^2)
%   p_Linf = ||p-ph||_inf,  p_L2=||p-ph||_0,  p_H1=|p-ph|_1
%
% L2、H1 的单元积分采用 7 点 5 阶三角高斯公式；L∞ 在单元内部细分采样求真上确界。
% [修改13] 原函数只计算弹性问题的向量位移误差（2*Nb 自由度）；现扩展为 Stokes
%         的“速度 + 压力”两套误差，压力段按一次元 201 在 Pb_p/Tb_p 上插值。

%% ---------- 自由度分段与局部基个数 ----------
Nb =size(Pb_u,2);
Nbp=size(Pb_p,2);
if numel(solution) ~= 2*Nb+Nbp
    error('compute_error_2D：solution 长度应为 2*Nb+Nbp（u1、u2、p）。');
end
sol1 =solution(1:Nb);            % u1 节点值
sol2 =solution(Nb+1:2*Nb);       % u2 节点值
solp =solution(2*Nb+1:end);      % p  节点值

switch basis_u
    case 201, nlb_u=3;
    case 202, nlb_u=6;
    otherwise, error('compute_error_2D：速度 basis_type 仅支持 201/202。');
end
switch basis_p
    case 201, nlb_p=3;
    otherwise, error('compute_error_2D：压力 basis_type 仅支持 201。');
end

Gpn=5;
number_of_elements=size(T,2);

%% ---------- 一、L2 与 H1 半范误差（高斯积分） ----------
L2_1=0; L2_2=0; H1_1=0; H1_2=0;   % 速度
L2_p=0; H1_p=0;                   % 压力

for n=1:number_of_elements

    vertices=P(:,T(:,n));
    [Gauss_weights,Gauss_nodes]=generate_Gauss_local_tri(vertices,Gpn);
    Ng=length(Gauss_weights);

    for k=1:Ng
        x=Gauss_nodes(k,1);
        y=Gauss_nodes(k,2);

        % ---- 速度数值解及其梯度（二次元局部基线性组合）----
        uh1=0; uh1x=0; uh1y=0;
        uh2=0; uh2x=0; uh2y=0;
        for beta=1:nlb_u
            phi  =FE_basis_local_fun_2D(x,y,vertices,basis_u,beta,0,0);
            phix =FE_basis_local_fun_2D(x,y,vertices,basis_u,beta,1,0);
            phiy =FE_basis_local_fun_2D(x,y,vertices,basis_u,beta,0,1);
            c1=sol1(Tb_u(beta,n));
            c2=sol2(Tb_u(beta,n));
            uh1 =uh1 +c1*phi;  uh1x=uh1x+c1*phix;  uh1y=uh1y+c1*phiy;
            uh2 =uh2 +c2*phi;  uh2x=uh2x+c2*phix;  uh2y=uh2y+c2*phiy;
        end

        % ---- 压力数值解及其梯度（一次元局部基线性组合）----
        ph=0; phx=0; phy=0;
        for beta=1:nlb_p
            phi  =FE_basis_local_fun_2D(x,y,vertices,basis_p,beta,0,0);
            phix =FE_basis_local_fun_2D(x,y,vertices,basis_p,beta,1,0);
            phiy =FE_basis_local_fun_2D(x,y,vertices,basis_p,beta,0,1);
            cp=solp(Tb_p(beta,n));
            ph =ph +cp*phi;  phx=phx+cp*phix;  phy=phy+cp*phiy;
        end

        % ---- 解析解及其梯度 ----
        u1 =exact_u_2D(x,y,1,0,0); u1x=exact_u_2D(x,y,1,1,0); u1y=exact_u_2D(x,y,1,0,1);
        u2 =exact_u_2D(x,y,2,0,0); u2x=exact_u_2D(x,y,2,1,0); u2y=exact_u_2D(x,y,2,0,1);
        pe =exact_u_2D(x,y,3,0,0); pex=exact_u_2D(x,y,3,1,0); pey=exact_u_2D(x,y,3,0,1);

        L2_1=L2_1+Gauss_weights(k)*(uh1-u1)^2;
        L2_2=L2_2+Gauss_weights(k)*(uh2-u2)^2;
        H1_1=H1_1+Gauss_weights(k)*((uh1x-u1x)^2+(uh1y-u1y)^2);
        H1_2=H1_2+Gauss_weights(k)*((uh2x-u2x)^2+(uh2y-u2y)^2);

        L2_p=L2_p+Gauss_weights(k)*(ph-pe)^2;
        H1_p=H1_p+Gauss_weights(k)*((phx-pex)^2+(phy-pey)^2);
    end
end

u_L2 =sqrt(L2_1+L2_2);
u_H1 =sqrt(H1_1+H1_2);
p_L2 =sqrt(L2_p);
p_H1 =sqrt(H1_p);

%% ---------- 二、L∞ 误差：Ω 上的真上确界（逐单元细分采样） ----------
n_sub=20;                                   % 每个单元每条参考方向细分数
[II,JJ]=meshgrid(0:n_sub,0:n_sub);
mask=(II+JJ)<=n_sub;                        % 只保留参考单元内的采样点
xi_s=(II(mask)/n_sub).';
eta_s=(JJ(mask)/n_sub).';

sup1=0; sup2=0; supp=0;
for n=1:number_of_elements

    vertices=P(:,T(:,n));
    % 参考坐标 -> 物理坐标（与 FE_basis_local_fun_2D 中的仿射映射一致）
    X=vertices(1,1)+(vertices(1,2)-vertices(1,1))*xi_s+(vertices(1,3)-vertices(1,1))*eta_s;
    Y=vertices(2,1)+(vertices(2,2)-vertices(2,1))*xi_s+(vertices(2,3)-vertices(2,1))*eta_s;

    uh1=0; uh2=0; ph=0;
    for beta=1:nlb_u
        phi=FE_basis_local_fun_2D(X,Y,vertices,basis_u,beta,0,0);
        uh1=uh1+sol1(Tb_u(beta,n))*phi;
        uh2=uh2+sol2(Tb_u(beta,n))*phi;
    end
    for beta=1:nlb_p
        phi=FE_basis_local_fun_2D(X,Y,vertices,basis_p,beta,0,0);
        ph=ph+solp(Tb_p(beta,n))*phi;
    end

    sup1 =max(sup1, max(abs(uh1-exact_u_2D(X,Y,1,0,0))));
    sup2 =max(sup2, max(abs(uh2-exact_u_2D(X,Y,2,0,0))));
    supp =max(supp, max(abs(ph -exact_u_2D(X,Y,3,0,0))));
end

u_Linf=max(sup1,sup2);
p_Linf=supp;

end
