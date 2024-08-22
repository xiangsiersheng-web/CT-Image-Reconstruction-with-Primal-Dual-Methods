% input：
%     noisy: 噪声图像
%     mu:正则化参数，一个常量
%     opts:一些参数设置，如步长、迭代次数等
% output：
%     noisyed：去噪后的图像
%     outs：一些输出参数
function [noisyed,outs]=hypd2_tv_rof(noisy,beita,im,opts)
    
    %% PDHG算法中的算子
    A = @(x) grad(x); % 向前差分的梯度算子
    At = @(y) -div(y); % 向后差分的负散度
    
    %% 迭代初始值
    [rows,cols] = size(noisy);
    x = noisy;
    y = A(noisy);
    %% opts是否存在、opts初始化
    if~exist('opts','var')
        opts = [];
    end
    opts = setDefaults(opts,noisy,A,At);
    
    %% 获取opts中的输入参数
    tau = opts.tau;        % 原问题步长
    sigma = opts.sigma;    % 对偶问题步长
    maxIters = opts.maxIters;
    theta = opts.theta;
    deltaT = opts.deltaT;
    gamma = opts.gamma;
                 
    
    %% 初始化输出参数outs
    outs = [];
    outs.L = opts.L; % A'A的谱半径的倒数
    outs.snr = zeros(maxIters+1,1);outs.snr(1) = 20*log10(norm(im,'fro')/norm(x-im,'fro'));
    outs.snrT = zeros(maxIters+1,2);

    t = 0; k = 1; outs.snrT(1,1)=0; outs.snrT(1,2) = outs.snr(1);
    %% 开始迭代
    for iter = 1:maxIters  
        tic;
        x0 = x;
        y0 = y;   
        %% 预测       
        % 对偶-更新
        ytilde = projectInf(y0+sigma*A(x0),1);
        ytilde = ytilde+theta*(ytilde-y0);
        % 原-更新
        xtilde = (beita*noisy+x0/tau-At(ytilde))/(beita+1/tau);
        
        %% 修正
        x = x0-gamma*(x0-xtilde);
        y = projectInf(y0-gamma*(y0-ytilde),1);        
        outs.snr(iter+1) = 20*log10(norm(im,'fro')/norm(x-im,'fro'));
        t=t+toc;
        if t>k*deltaT
            outs.snrT(k+1,1)=k*deltaT; outs.snrT(k+1,2) = outs.snr(iter+1);
            k = k+1;
        end
        % 迭代停止条件
%         if iter >= maxIters || (primal/maxPrimal<opts.tol && dual/maxDual<opts.tol && iter>5) || (primal<1e-10 && dual<1e-10)
        if iter >= maxIters || t>opts.T
            outs.snr = outs.snr(1:iter+1);
            outs.snrT = outs.snrT(1:k,:);
            outs.iters = iter;
            noisyed = x;
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

function z = projectInf(z,lambda)
x = z(:,:,1);
y = z(:,:,2);
norm = max(sqrt(x.*x+y.*y),lambda);
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

function opts = setDefaults(opts,x0,A,At)
%  L:  A'A谱半径的倒数，sigma*tau不大于L，否则不能收敛
if ~isfield(opts,'L') || opts.L<=0
%     x = randn(size(x0));
%     transform = At(A(x));
%     specRadius = norm(transform,'fro')/norm(x,'fro');
%     opts.L = 2/specRadius;
    s = NSV(A,At,x0,1e-8,1e4);
    opts.L = 1/s;
end

%  maxIters: 最大迭代次数
if ~isfield(opts,'maxIters')
    opts.maxIters = 1e4;
end
% tol:  残差的相对减小<tol时，停止
if ~isfield(opts,'tol') 
    opts.tol = 1e-3;
end

if ~isfield(opts,'theta')        
    opts.theta = 1;
end

if ~isfield(opts,'T')        
    opts.T = 100;
end

if ~isfield(opts,'deltaT')        
    opts.deltaT = 0.1;
end

if ~isfield(opts,'gamma')        
    opts.gamma = 1.9;
end

% tau:  原问题迭代步长
if ~isfield(opts,'tau')        
    opts.tau = sqrt(opts.L);
end
% sigma: 对偶问题的迭代步长
if ~isfield(opts,'sigma')      
    opts.sigma = opts.L/opts.tau;
end
return