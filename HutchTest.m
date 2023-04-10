rng(1);
n   = 128;
%A   = randn(n,n);
A = logm(M);
% e.g., A = log(M)

nv  = 1e5;
V   = sign(randn(n,nv));
AV  = A*V;
VtAV= sum( V .* AV , 1);

er = abs(cumsum( VtAV )./(1:nv) - trace(A));
figure(1)
loglog( 1:nv, er )