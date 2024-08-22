%% REFERENCE
% https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

%% ART Equation
% x^(k+1) = x^k + lambda * AT(b - A(x))/ATA

%%
function [x,outs]  = ART_test(A,AT,b,im,lambda,MaxIt)
if (nargin < 6)
    MaxIt       = 1e3;
end

x = zeros(size(im));
outs.snr = zeros(MaxIt+1,1);outs.snr(1) = 20*log10(norm(im,'fro')/norm(x-im,'fro'));
outs.snrT = zeros(MaxIt+1,2);outs.snrT(1,1) = 0;outs.snrT(1,2) = outs.snr(1);
outs.toll = zeros(MaxIt,1);
outs.tol = zeros(MaxIt+1,1);outs.tol(1) = norm(A(x)-b,'fro')/length(x(:));
outs.tol111 = zeros(MaxIt+1,1);outs.tol111(1) = ((1/2)*norm(A(x)-b,'fro')^2+TV(grad(x)))/length(x(:));

ATA	= AT(A(ones(size(x), 'single')));

deltaT=0.1; t = 0; Tk = 1;
for k = 1:MaxIt
    tic;
    x0 = x;
    x	= x0 + lambda*AT(b - A(x0))./ATA;  
    x(x < 0) = 0;
    outs.snr(k+1) = 20*log10(norm(im,'fro')/norm(x-im,'fro')); 
    t = t+toc;
    if t>Tk*deltaT
        outs.snrT(Tk+1,1) = Tk*deltaT;
        outs.snrT(Tk+1,2) = outs.snr(k+1);
        Tk = Tk+1;
    end 
    outs.tol(k+1) = norm(A(x)-b,'fro')/length(b(:));
    outs.tol111(k+1) = ((1/2)*norm(A(x)-b,'fro')^2+TV(grad(x)))/length(x(:));
    outs.toll(k) = norm(x-x0,'fro')/length(x(:));
    if k>=MaxIt
        outs.tol = outs.tol(1:k+1);
        outs.tol111 = outs.tol111(1:k+1);
        outs.toll = outs.toll(1:k);
        outs.snr = outs.snr(1:k+1);
        outs.snrT = outs.snrT(1:Tk,:);
        return;
    end
%     figure(1); colormap gray;
%     imagesc(x);
%     axis image off;
%     title(['迭代次数',num2str([i, niter], '%d / %d')]);
%     drawnow();
end
return

function g = grad(u)
 [rows,cols] = size(u);
 g = zeros(rows,cols,2);
 g(:,:,1) = Dx(u);
 g(:,:,2) = Dy(u);
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