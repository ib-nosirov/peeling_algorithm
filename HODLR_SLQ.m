% sample matrix
n=400; diagSize=125; k=10; I=[1 n];
A = orth(randn(n)); A = A*diag(n+1:2*n)*A';
A = MakeHODLRMtrx(A,n,k,diagSize,I);
%A = randn(n);
%A = A * A';
q = randn(n,1);
m = 20;
nv = 10;
Gamma = 0;
for ii=1:nv
    v = 2*round(rand(n,1))-1; v = v/norm(v);
    T = HODLR_Lanczos(A,v,m+1);
    [Y,Theta] = eig(T);
    t = Y(1,:);
    t = t.^2;
    theta = diag(Theta);
    theta = log(theta);
    Gamma = Gamma + sum(t.*theta');
end
Gamma = n/nv * Gamma