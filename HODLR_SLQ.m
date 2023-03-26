function Gamma = HODLR_SLQ(A,n,f,m,nv)
q = randn(n,1); % check the code
Gamma = 0;
for ii=1:nv
    v = 2*round(rand(n,1))-1; v = v/norm(v);
    T = HODLR_Lanczos(A,v,m);
    [Y,Theta] = eig(T);
    theta = diag(Theta);
    tau = Y(1,:);
    % Idea, remove the corresponding entries. (play with the threshold,
    % should be very small)
    tau = tau(theta > 1e-10).^2;
    theta = f(theta(theta > 1e-10));

    Gamma = Gamma + sum(tau.*theta');
end
Gamma = n/nv * Gamma;
