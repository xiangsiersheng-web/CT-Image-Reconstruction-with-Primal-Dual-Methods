% input：
%     f:redon变化后的数据
%     W:radon变换矩阵
%     beita:正则化参数，一个常量
%     rows,cols:图像的行列数
%     im:用于计算信噪比
%     opts:一些参数设置，如步长、迭代次数等
% output：
%     reimed：去噪后的图像
%     outs：一些输出参数
function [reimed, outs] = cp_tv_me(f,W,beita,rows,cols,im,opts)
    %% 迭代初始值
    y = zeros(rows,cols);
    x = zeros(size(f));
    z = zeros(rows,cols,2);
    yy = zeros(rows,cols);
    
    %% CP算法中的算子
    A = @(y) grad(y); % 向前差分的梯度算子
    At = @(z) -div(z); % 向后差分的负散度
    
    %% opts是否存在、opts初始化
    if~exist('opts','var')
        opts = [];
    end
    opts = setDefaults(opts,y,A,At,W);
    
    %% 获取opts中的输入参数
    tau = opts.tau;    % 对偶问题步长
    sigma = opts.sigma;        % 原问题步长
    maxIters = opts.maxIters;
    tol = opts.tol;
    theta = opts.theta;
    deltaT = opts.deltaT;
    
    %% 初始化输出参数outs
    outs = [];
    outs.L = opts.L;             % K'K的谱半径的倒数
    outs.snr = zeros(maxIters+1,1);
    outs.tol = zeros(maxIters+1,1);
    outs.snr(1) = 20*log10(norm(im,'fro')/norm(y-im,'fro'));
%     outs.tol(1) = norm(W*y(:)-f)/length(f);
%     outs.tol(1) = (beita/2)*norm(W*y(:)-f)^2/length(f)+TV(A(y))/length(y(:));
    outs.tol(1) = ((beita/2)*norm(W*y(:)-f)^2+TV(A(y)))/length(y(:));
    outs.snrT = zeros(maxIters,2);
    
    t=0;k=1;%用于记录时间，时间间隔
    outs.snrT(1,1) = t;outs.snrT(1,2)=outs.snr(1);  
    %% 开始迭代
    for iter = 1:maxIters   
        tic;
        %% 保存旧值
        y0 = y;
        x0 = x;
        z0 = z;
        yy0 = yy;
                
        %% 对偶-更新
        z = projectInf(z0+tau*A(yy0),1);
        x = (x0+tau*(W*yy0(:)-f))/(1+tau/beita);

        %% 原-更新
        y = y0 - sigma*(reshape(W'*x,rows,cols) + At(z));
        yy = y + theta*(y - y0);

        
        
        
        %% 按照时间来取信噪比        
        t=t+toc;
        %% 计算信噪比
        outs.snr(iter+1) = 20*log10(norm(im,'fro')/norm(y-im,'fro'));
        if  t>k*deltaT
            outs.snrT(k+1,1)=k*deltaT;
            outs.snrT(k+1,2)=outs.snr(iter+1);
            k = k+1;
        end
%         outs.tol(iter+1) = norm(W*y(:)-f)/length(f);
%         outs.tol(iter+1) = (beita/2)*norm(W*y(:)-f)^2/length(f)+TV(A(y))/length(y(:));
        outs.tol(iter+1) = ((beita/2)*norm(W*y(:)-f)^2+TV(A(y)))/length(y(:));
        %% 迭代停止条件
        if  outs.snr(iter+1)>opts.snr||iter>=maxIters||outs.tol(iter+1)<tol||outs.snrT(k,1)>opts.T
            outs.iters = iter;
            outs.snr = outs.snr(1:iter+1);
            outs.snrT = outs.snrT(1:k,:);
            outs.tol = outs.tol(1:iter+1);
            reimed = y;           
            return;
        end
        
    end
return

function g = grad(u)
 [rows,cols] = size(u);
 g = zeros(rows,cols,2);
 g(:,:,1) = Dx(u);
 g(:,:,2) = Dy(u);
return
 
function d = div(u)
 d = -(Dxt(u(:,:,1))+Dyt(u(:,:,2)));
return

function z = projectInf(z,lamda)
x = z(:,:,1);
y = z(:,:,2);
norm = max(sqrt(x.*x+y.*y),lamda);
% norm = sqrt(x.*x+y.*y);
z(:,:,1) = x./norm;
z(:,:,2) = y./norm;
return

function d = Dx(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
return

function d = Dy(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
return

function d = Dxt(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
d(:,cols) = u(:,cols)-u(:,1);
return

function d = Dyt(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);
d(rows,:) = u(rows,:)-u(1,:);
return

function tv = TV(x)
[rows,cols,~]=size(x);
tv = 0;
for i = 1:rows
    for j = 1:cols
        c = [x(i,j,1),x(i,j,2)];
        tv = tv+norm(c);
    end
end
return

function opts = setDefaults(opts,y0,A,At,W)
%  L:  K'K谱半径的倒数，tau*sigma不大于L，否则不能收敛
if ~isfield(opts,'L') || opts.L<=0
    s = NSV(A,At,W,y0,1e-10,1e4);
    opts.L = 1/s;
end

%  deltaT: 时间间隔，默认0.1秒取一次信噪比
if ~isfield(opts,'deltaT')
    opts.deltaT = 0.1;
end

%  T: 运行时间限制，默认1000秒
if ~isfield(opts,'T')
    opts.T = 1000;
end

% theta: 算法中一个参数
if ~isfield(opts,'theta')
    opts.theta = 1;
end
%  maxIters: 最大迭代次数
if ~isfield(opts,'maxIters')
    opts.maxIters = 1e4;
end

% tol:  残差的相对减小<tol时，停止
% if ~isfield(opts,'tol') 
%     opts.tol = 1e-3;
% end

% tol:迭代停止
if ~isfield(opts,'tol')
    opts.tol = 1e-6;
end

% snr:迭代停止
if ~isfield(opts,'snr')
    opts.snr = 30;
end


if ~isfield(opts,'miu')
    opts.miu = 1;
end
% sigma:  原问题迭代步长
if ~isfield(opts,'sigma')        
    opts.sigma = sqrt(opts.miu*opts.L);
end
% tau: 对偶问题的迭代步长
if ~isfield(opts,'tau')      
    opts.tau = opts.L*opts.miu/opts.sigma;
end
return

%幂法求K的奇异值
function s = NSV(A,At,W,y0,tol,N)
[rows,cols] = size(y0);
y0 = randn(rows,cols);
wtwy = reshape(W'*(W*y0(:)),rows,cols);
y1 = At(A(y0))+wtwy;
n = 1;
while abs((max(max(y0))-max(max(y1))))>tol && n<=N
    y0 = y1;
    wtwy = reshape(W'*(W*y0(:)/max(max(y0))),rows,cols);
    y1 = At(A(y0/max(max(y0))))+wtwy;
    n = n+1;
end
s = max(max(y1));
return;
