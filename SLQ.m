% pass in function handle for A, rather than A.
function Gamma = SLQ(kMtrxFcn,n,f,m,nv)
q = randn(n,1);
Gamma = 0;
for ii=1:nv
    v = 2*round(rand(n,1))-1; v = v/norm(v);
    T = Lanczos(kMtrxFcn,v,m+1);
    [Y,Theta] = eig(T);
    t = Y(1,:);
    t = sort(t.^2);

    theta = sort(diag(Theta));
    theta = f(theta(theta > 0.1));
    Gamma = Gamma + sum(t(length(t)+1-length(theta):end).*theta');
end
Gamma = n/nv * Gamma;
