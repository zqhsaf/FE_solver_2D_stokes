function [Pb,Tb]=generate_PbTb(left,right,bottom,top,N1,N2,basis_type)

if basis_type==201
    [Pb,Tb]=generate_PT(left,right,bottom,top,N1,N2);

elseif basis_type==202
    h1=(right-left)/N1;
    h2=(top-bottom)/N2;

    M=2*N1+1;          % x 方向二次节点数
    L=2*N2+1;          % y 方向二次节点数
    Nq=M*L;            % 二次节点总数
    Ne=2*N1*N2;        % 三角单元总数

    %% 1. 统一规则网格节点坐标（按列优先编号）
    Pb=zeros(2,Nq);
    for I=0:2*N1
        for J=0:2*N2
            k=I*L+J+1;
            Pb(1,k)=left+I*h1/2;
            Pb(2,k)=bottom+J*h2/2;
        end
    end

    % 二次网格 (I,J) -> 全局编号
    qid=@(I,J) I*L+J+1;

    %% 2. 三角单元连接表（与 generate_PT 的单元顺序、朝向一致）
    Tb=zeros(6,Ne);
    for i=0:N1-1
        for j=0:N2-1
            n=2*(i*N2+j)+1;

            % 左下三角：顶点 (i,j)、(i+1,j)、(i,j+1)
            %   局部4=边1-2（水平）中点，局部5=边2-3（对角线）中点，局部6=边1-3（竖直）中点
            Tb(:,n)=[ ...
                qid(2*i,   2*j);
                qid(2*i+2, 2*j);
                qid(2*i,   2*j+2);
                qid(2*i+1, 2*j);
                qid(2*i+1, 2*j+1);
                qid(2*i,   2*j+1)];

            % 右上三角：顶点 (i,j+1)、(i+1,j)、(i+1,j+1)
            %   局部4=边1-2（对角线）中点，局部5=边2-3（竖直）中点，局部6=边1-3（水平）中点
            Tb(:,n+1)=[ ...
                qid(2*i,   2*j+2);
                qid(2*i+2, 2*j);
                qid(2*i+2, 2*j+2);
                qid(2*i+1, 2*j+1);
                qid(2*i+2, 2*j+1);
                qid(2*i+1, 2*j+2)];
        end
    end

else
    error('generate_PbTb：不支持的 basis_type（201 一次，202 二次）。');
end

end