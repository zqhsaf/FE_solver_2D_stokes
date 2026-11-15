function run_error_table
% 复现讲义第 52、53 页的两张收敛性误差表（Stokes Example 1）
%   区域 Ω=[0,1]×[-0.25,0]；速度二次元 202，压力一次元 201（Taylor-Hood）
%   表 1（第 52 页）：速度误差 ||u-uh||_inf、||u-uh||_0、|u-uh|_1
%   表 2（第 53 页）：压力误差 ||p-ph||_inf、||p-ph||_0、|p-ph|_1
% 网格步长 h = 1/8、1/16、1/32、1/64（N1=N2=8,16,32,64）。
% [修改14] 原脚本调用不存在的 FE_solver_2D_elasticity、区域为 [0,1]×[0,1]，
%         且只输出位移误差；现改为 Stokes 求解器，并同时输出速度、压力两张表。

clc;
left=0; right=1; bottom=-0.25; top=0;
N_all =[8 16 32 64];

basis_u=202;   % 速度：二次三角元
basis_p=201;   % 压力：一次三角元

Eu=zeros(numel(N_all),3);   % 速度：Linf, L2, H1
Ep=zeros(numel(N_all),3);   % 压力：Linf, L2, H1

for k=1:numel(N_all)
    N=N_all(k);
    t0=tic;
    [sol,P,T,Pb_u,Tb_u,Pb_p,Tb_p]=FE_solver_2D_stokes(left,right,bottom,top,N,N);
    [uLinf,uL2,uH1,pLinf,pL2,pH1]=compute_error_2D(sol,P,T, ...
        Pb_u,Tb_u,basis_u,Pb_p,Tb_p,basis_p);
    Eu(k,:)=[uLinf,uL2,uH1];
    Ep(k,:)=[pLinf,pL2,pH1];
    fprintf('h=1/%-3d  求解+误差用时 %.1f s\n',N,toc(t0));
end

%% ---- 表 1：速度误差（讲义第 52 页）----
fprintf('\n===== 表 1：速度（二次三角元 202）数值误差（讲义第 52 页）=====\n');
fprintf('%8s  %14s  %14s  %14s\n','h','||u-uh||_inf','||u-uh||_0','|u-uh|_1');
for k=1:numel(N_all)
    fprintf('%8s  %14.4e  %14.4e  %14.4e\n', ...
        sprintf('1/%d',N_all(k)),Eu(k,1),Eu(k,2),Eu(k,3));
end
fprintf('速度收敛阶 log2(e_{k-1}/e_k)：\n');
for k=2:numel(N_all)
    fprintf('  h=1/%-3d -> 1/%-3d :  %.2f   %.2f   %.2f\n', ...
        N_all(k-1),N_all(k), ...
        log2(Eu(k-1,1)/Eu(k,1)),log2(Eu(k-1,2)/Eu(k,2)),log2(Eu(k-1,3)/Eu(k,3)));
end

%% ---- 表 2：压力误差（讲义第 53 页）----
fprintf('\n===== 表 2：压力（一次三角元 201）数值误差（讲义第 53 页）=====\n');
fprintf('%8s  %14s  %14s  %14s\n','h','||p-ph||_inf','||p-ph||_0','|p-ph|_1');
for k=1:numel(N_all)
    fprintf('%8s  %14.4e  %14.4e  %14.4e\n', ...
        sprintf('1/%d',N_all(k)),Ep(k,1),Ep(k,2),Ep(k,3));
end
fprintf('压力收敛阶 log2(e_{k-1}/e_k)：\n');
for k=2:numel(N_all)
    fprintf('  h=1/%-3d -> 1/%-3d :  %.2f   %.2f   %.2f\n', ...
        N_all(k-1),N_all(k), ...
        log2(Ep(k-1,1)/Ep(k,1)),log2(Ep(k-1,2)/Ep(k,2)),log2(Ep(k-1,3)/Ep(k,3)));
end

%% ---- 结果写入文本文件（防覆盖：带时间戳）----
ts=datestr(now,'yyyymmdd_HHMMSS');
outfile=fullfile(fileparts(mfilename('fullpath')),['收敛性表格运行结果_',ts,'.txt']);
fid=fopen(outfile,'w','n','UTF-8');
fprintf(fid,'Stokes Taylor-Hood 有限元收敛性误差表（讲义第 52、53 页）\n');
fprintf(fid,'区域 Ω=[0,1]×[-0.25,0]，速度 202 / 压力 201，ν=1\n\n');
fprintf(fid,'表 1：速度误差（第 52 页）\n');
fprintf(fid,'%8s  %14s  %14s  %14s\n','h','||u-uh||_inf','||u-uh||_0','|u-uh|_1');
for k=1:numel(N_all)
    fprintf(fid,'%8s  %14.4e  %14.4e  %14.4e\n', ...
        sprintf('1/%d',N_all(k)),Eu(k,1),Eu(k,2),Eu(k,3));
end
fprintf(fid,'\n表 2：压力误差（第 53 页）\n');
fprintf(fid,'%8s  %14s  %14s  %14s\n','h','||p-ph||_inf','||p-ph||_0','|p-ph|_1');
for k=1:numel(N_all)
    fprintf(fid,'%8s  %14.4e  %14.4e  %14.4e\n', ...
        sprintf('1/%d',N_all(k)),Ep(k,1),Ep(k,2),Ep(k,3));
end
fclose(fid);
fprintf('\n结果已保存至：%s\n',outfile);

end
