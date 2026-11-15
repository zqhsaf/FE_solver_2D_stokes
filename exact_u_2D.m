function value=exact_u_2D(x,y,comp,der_x,der_y)
% 二维 Stokes 算例（讲义第 50 页）的解析解及其偏导数
%   comp=1：u1 = x^2 y^2 + e^{-y}
%   comp=2：u2 = -(2/3) x y^3 + 2 - pi sin(pi x)
%   comp=3： p = -(2 - pi sin(pi x)) cos(2 pi y)
% der_x,der_y：对 x、y 的偏导阶数（只支持 0/1 阶，更高阶返回 0）
% 支持 x、y 为同维数组（逐元素运算），便于在节点上向量化求值。
%
% [修改2] 原文件是二维弹性方程解析解（sin(pi x)sin(pi y) 等），与本 Stokes
%         算例不符；现按讲义第 50 页解析解整体替换，并新增 comp=3 的压力 p。

if comp==1
    % ---- u1 = x^2 y^2 + e^{-y} ----
    if der_x==0 && der_y==0
        value=x.^2.*y.^2 + exp(-y);
    elseif der_x==1 && der_y==0
        value=2*x.*y.^2;                       % ∂u1/∂x
    elseif der_x==0 && der_y==1
        value=2*x.^2.*y - exp(-y);             % ∂u1/∂y
    else
        value=zeros(size(x));
    end

elseif comp==2
    % ---- u2 = -(2/3) x y^3 + 2 - pi sin(pi x) ----
    if der_x==0 && der_y==0
        value=-(2/3)*x.*y.^3 + 2 - pi*sin(pi*x);
    elseif der_x==1 && der_y==0
        value=-(2/3)*y.^3 - pi^2*cos(pi*x);    % ∂u2/∂x
    elseif der_x==0 && der_y==1
        value=-2*x.*y.^2;                       % ∂u2/∂y
    else
        value=zeros(size(x));
    end

elseif comp==3
    % ---- p = -(2 - pi sin(pi x)) cos(2 pi y) ----
    if der_x==0 && der_y==0
        value=-(2 - pi*sin(pi*x)).*cos(2*pi*y);
    elseif der_x==1 && der_y==0
        value=pi^2*cos(pi*x).*cos(2*pi*y);     % ∂p/∂x
    elseif der_x==0 && der_y==1
        value=2*pi*(2 - pi*sin(pi*x)).*sin(2*pi*y);  % ∂p/∂y
    else
        value=zeros(size(x));
    end

else
    error('exact_u_2D：comp 只能取 1（u1）、2（u2）或 3（p）。');
end

end
