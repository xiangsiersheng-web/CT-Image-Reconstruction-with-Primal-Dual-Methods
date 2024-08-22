%% REFERENCE
% https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

%% ART Equation
% x^(k+1) = x^k + lambda * AT(b - A(x))/ATA

%%
function x  = art(A,b,lambda,niter)

[M,N] = size(A);
x = ones(N,1);

% ATA	= A'*(A*ones(N,1));
% % t1 = 0;t2 = 0;
% for k = 1:niter
% %     tic;
%     x0 = x;
%     bPj = A'*(b - A*x0);
% %     t1=t1+toc;tic;
%     for i = 1:N
%         x(i) = x0(i) + lambda*bPj(i)./ATA(i);
%         if x(i)<0
%             x(i) = 0;
%         end
%     end
% %     figure(1); colormap gray;
% %     utem = reshape(x,sqrt(N),sqrt(N));
% %     imagesc(utem);
% %     axis image off;
% %     title(['迭代次数:',num2str([k, niter], '%d / %d')]);
% %     drawnow();
% %     t2 =t2+toc;
% end
% t1
% t2

norm2_row = zeros(M,1);
Pos = zeros(M,N);
for j = 1:M
    norm2_row(j) = norm(A(j,:))^2;
    tem = find(A(j,:)>0);
    Pos(j,1:length(tem))=tem;
end
t1=0;
for k = 1:niter
    x0 = x;
    bPj = (b-A*x0).*A;
    for i = 1:M 
        for j = nonzeros(Pos(i,:))'
            if isempty(j)
                continue;
            end
            tt1 =tic;
            x(j) = x0(j) + lambda*bPj(i,j)/norm2_row(i);             
            if x(j)<0
                x(j) = 0;
            end
            t1=t1+toc(tt1);
        end
    end      
end
t1
% t1 = 0;t2 =0;
% for k = 1:niter
%     x0 = x;
%     bPj = (b-A*x0).*A;
%     for i = 1:M 
%         tt1 = tic;
%         for j = 1:N
%             if bPj(i,j)~=0
%                 tt2 = tic;
%                 x(j) = x(j) + lambda*bPj(i,j)/norm2_row(i); 
%                 if x(j)<0
%                     x(j) = 0;
%                 end
%                 t2 = t2+toc(tt2);
%             end
%         end
%         t1 = t1+toc(tt1);
%     end  
% end
% t1
% t2
end