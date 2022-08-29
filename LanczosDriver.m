% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
n = 1000; diagSize= 125; k = 10; I = [1 n^d];
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
A = rbf(ep,DM);
q = randn(n,1);
m = 20;
T = Lanczos(A,q,m)
[V,D] = eig(T);
%    sort eigenvalues into decreasing order
 [Ds,Is] = sort(-diag(D));
 Eval = -Ds;
 Eval