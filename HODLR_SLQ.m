function Gamma = HODLR_SLQ(A,n,f,m,nv)
q = randn(n,1);
Gamma = 0;
for ii=1:nv
    v = 2*round(rand(n,1))-1; v = v/norm(v);
    T = HODLR_Lanczos(A,v,m);
    [Y,Theta] = eig(T);
    t = Y(1,:);
    t = sort(t.^2);

    theta = sort(diag(Theta));
    theta = f(theta(theta > 0.1));
    Gamma = Gamma + sum(t(length(t)+1-length(theta):end).*theta');
end
Gamma = n/nv * Gamma;
