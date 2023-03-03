% pass in function handle for A, rather than A.
function Gamma = SLQ(kMtrxFcn,n,f,m,nv)
q = randn(n,1);
Gamma = 0;
for ii=1:nv
    v = 2*round(rand(n,1))-1; v = v/norm(v);
    T = Lanczos(kMtrxFcn,v,m);
    [Y,Theta] = eig(T);
    theta = diag(Theta);
    tau = Y(1,:);
    % Idea, remove the corresponding entries.
    tau = tau(theta > 0.1).^2;
    theta = f(theta(theta > 0.1));

    Gamma = Gamma + sum(tau.*theta');
end
Gamma = n/nv * Gamma;
