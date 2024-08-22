%幂法求A的奇异值
function s = NSV(A,At,x0,tol,N)
x0 = randn(size(x0));
x1 = At(A(x0));
n = 1;
while abs((max(max(x0))-max(max(x1))))>tol && n<=N
    x0 = x1;
    x1 = At(A(x0./max(max(x0))));
    n = n+1;
end
s = max(max(x1));
return;