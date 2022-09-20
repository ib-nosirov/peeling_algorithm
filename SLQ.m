function Gamma = SLQ(A,n,f,m,nv)
q = randn(n,1);
Gamma = 0;
for ii=1:nv
    v = 2*round(rand(n,1))-1; v = v/norm(v);
    T = Lanczos(A,v,m+1);
    [Y,Theta] = eig(T);
    t = Y(1,:);
    t = t.^2;
    theta = diag(Theta);
    theta = f(theta);
    Gamma = Gamma + sum(t.*theta');
end
Gamma = n/nv * Gamma;
