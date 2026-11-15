function plot_stokes(N)
% 绘制 Stokes 算例的数值解与解析解对比云图（速度大小、压力）
% 用法：plot_stokes(16)  （N 为 x/y 方向剖分数，默认 16）
% 图片高清保存到项目目录，文件名带时间戳，防止覆盖。
% [修改15] 新增可视化脚本：速度大小 |u|、压力 p 的数值解/解析解 2×2 对比。
% [修改16] 修复：速度为二次元节点、压力为一次元节点，二者连接表不同；不能用粗网格
%         T 直接 trisurf 二次元节点。统一改用 scatteredInterpolant 插值到规则网格，
%         再用 contourf 填色，避免节点-单元对应错乱导致的图形错位。

if nargin<1 || isempty(N), N=16; end
left=0; right=1; bottom=-0.25; top=0;

[solution,P,T,Pb_u,Tb_u,Pb_p,Tb_p]=FE_solver_2D_stokes(left,right,bottom,top,N,N);
Nb =size(Pb_u,2);  Nbp=size(Pb_p,2);
u1h=solution(1:Nb);            u2h=solution(Nb+1:2*Nb);   ph=solution(2*Nb+1:end);

% 解析解在节点上的值
u1e=exact_u_2D(Pb_u(1,:),Pb_u(2,:),1,0,0)';
u2e=exact_u_2D(Pb_u(1,:),Pb_u(2,:),2,0,0)';
pe =exact_u_2D(Pb_p(1,:),Pb_p(2,:),3,0,0)';

% 规则绘图网格
[Xg,Yg]=meshgrid(linspace(left,right,120),linspace(bottom,top,60));

% 数值解插值到规则网格（坐标与值统一为列向量）
F1h=scatteredInterpolant(Pb_u(1,:)',Pb_u(2,:)',u1h);
F2h=scatteredInterpolant(Pb_u(1,:)',Pb_u(2,:)',u2h);
Fph=scatteredInterpolant(Pb_p(1,:)',Pb_p(2,:)',ph);
u1hg=F1h(Xg,Yg); u2hg=F2h(Xg,Yg); phg=Fph(Xg,Yg);
mag_h=sqrt(u1hg.^2+u2hg.^2);

% 解析解在规则网格上直接求值
u1eg=exact_u_2D(Xg,Yg,1,0,0);
u2eg=exact_u_2D(Xg,Yg,2,0,0);
peg =exact_u_2D(Xg,Yg,3,0,0);
mag_e=sqrt(u1eg.^2+u2eg.^2);

set(0,'DefaultAxesFontName','Microsoft YaHei');
figure('Color','w','Position',[100 100 1100 800]);

subplot(2,2,1); contourf(Xg,Yg,mag_h,20,'LineStyle','none'); colorbar;
title('数值解：速度大小 |u_h|'); xlabel('x'); ylabel('y'); axis equal tight;

subplot(2,2,2); contourf(Xg,Yg,mag_e,20,'LineStyle','none'); colorbar;
title('解析解：速度大小 |u|'); xlabel('x'); ylabel('y'); axis equal tight;

subplot(2,2,3); contourf(Xg,Yg,phg,20,'LineStyle','none'); colorbar;
title('数值解：压力 p_h'); xlabel('x'); ylabel('y'); axis equal tight;

subplot(2,2,4); contourf(Xg,Yg,peg,20,'LineStyle','none'); colorbar;
title('解析解：压力 p'); xlabel('x'); ylabel('y'); axis equal tight;

sgtitle(sprintf('Stokes Taylor-Hood 有限元解对比（h=1/%d）',N));

ts=datestr(now,'yyyymmdd_HHMMSS');
outfile=fullfile(fileparts(mfilename('fullpath')),['Stokes解对比图_h',num2str(N),'_',ts,'.png']);
exportgraphics(gcf,outfile,'Resolution',300);
fprintf('图像已保存：%s\n',outfile);

end
