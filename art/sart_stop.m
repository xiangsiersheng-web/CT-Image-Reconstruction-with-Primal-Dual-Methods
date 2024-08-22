%% REFERENCE
% https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

%%
function [v,outs] = sart_stop(b,W,lambda,MaxIt,im)
N = size(W,2);
v = zeros(N,1);
outs.snr = zeros(MaxIt+1,1);outs.snr(1) = 20*log10(norm(im,'fro')/norm(v-im(:)));
outs.snrT = zeros(MaxIt+1,2);outs.snrT(1,1) = 0;outs.snrT(1,2) = outs.snr(1);
outs.toll = zeros(MaxIt,1);
outs.tol = zeros(MaxIt+1,1);outs.tol(1) = norm(W*v-b)/length(b);
sumCol = full(sum(W)');       
sumRow = full(sum(W,2));
deltaT=0.1; t = 0; Tk = 1;
for k = 1:MaxIt
    tic;
    v0 = v;    
    bPj = W'* ((b-W*v0)./sumRow);
    for i=1:N
        v(i) = v0(i)+ lambda*bPj(i)/sumCol(i);
        if v(i)<0
            v(i) = 0;
        end
    end
    outs.snr(k+1) = 20*log10(norm(im,'fro')/norm(v-im(:)));  
    t = t+toc;      
    if t>Tk*deltaT
        outs.snrT(Tk+1,1) = Tk*deltaT;
        outs.snrT(Tk+1,2) = outs.snr(k+1);
        Tk = Tk+1;
    end
    outs.tol(k+1) = norm(W*v-b)/length(b);
    outs.toll(k) = norm(v-v0)/length(v);
    if k>=MaxIt||outs.snr(k+1)<outs.snr(k)
        outs.tol = outs.tol(1:k+1);
        outs.toll = outs.toll(1:k);
        outs.snr = outs.snr(1:k+1);
        outs.snrT = outs.snrT(1:Tk,:);
        return;
    end
end
% t1
% t2

end
