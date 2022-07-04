% make a kernel
N = 128;
K = zeros(N);
sigma = N/10;
for ii = 1:N
    for jj = 1:N
        K(ii,jj) = exp(-(ii-jj)^2/sigma^2);
    end
end
n = 128; diag_size = 8; k = 5; I = [1 n^d];
%rbf = @(e,r) exp(-r.^2/ep^2); ep = n/10; d = 1;
% Point evals at which to sample
%x = CreatePoints(n^d,d,'u');
% compute the absolute difference
%DM = DistanceMatrix(x,x);
% sample matrix
%K_mtrx = rbf(ep,DM);
%imagesc(K_mtrx);
exact = trace(K)
[U_tree Z_tree idx_tree] = NumericalExperiment(K,n^d,k,diag_size,I);

K_trace = HutchPlusPlus(K_mtrx,U_tree,idx_tree)
