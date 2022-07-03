% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 19; d = 1;
n = 1000; diag_size = 130; k = 5; I = [1 n^d];
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
K_mtrx = rbf(ep,DM);
%imagesc(K_mtrx);
exact = trace(K_mtrx);
[U_tree Z_tree idx_tree] = NumericalExperiment(K_mtrx,n^d,k,diag_size,I);

K_trace = HutchPlusPlus(K_mtrx,U_tree,idx_tree)
