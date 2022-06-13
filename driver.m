n = 100; diag_size = 13; k = 10; I = [1 n];
% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
K_mtrx = rbf(ep,DM);
[U_tree Z_tree idx_tree] = NumericalExperiment(K_mtrx,n,k,diag_size,I);
