% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
n = 1000; diagSize= 125; k = 10; I = [1 n^d];
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
n=400; A = orth(randn(n)); A = A*diag([n+1:2*n])*A'; 
%A = randn(n);
%A = A * A';
q = randn(n,1);
m = 20;
nv = 10;
Gamma = 0;
for ii=1:nv
    v = 2*round(rand(n,1))-1; v = v/norm(v);
    T = Lanczos(A,v,m+1);
    [Y,Theta] = eig(T);
    t = Y(1,:);
    t = t.^2;
    theta = diag(Theta);
    theta = log(theta);
    Gamma = Gamma + sum(t.*theta');
end
Gamma = n/nv * Gamma
