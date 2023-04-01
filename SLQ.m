% pass in function handle for A, rather than A.
function Gamma = SLQ(linearOperator,n,f,m,nv)
Gamma = 0;
for ii=1:nv
    % make a Rademacher random variable.
    v = 2*round(rand(n,1))-1;
    v = v/norm(v);
    T = Lanczos(linearOperator,v,m);
    [Y,Theta] = eig(T);
    theta = diag(Theta);
    tau = Y(1,:);
    % Idea, remove the corresponding entries.
    tau = tau(theta > 0.000001).^2;
    theta = f(theta(theta > 0.000001));
    Gamma = Gamma + sum(tau.*theta');
end
Gamma = n/nv*Gamma;
