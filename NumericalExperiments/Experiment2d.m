% Goal: test the accuracy of HODLR-SLQ in recovering a 2D kernel with
% regularized eigenvalues.

% Create a 2d C6 Matern kernel
kernel = @(e,r) (1+e*r+2/5*(e*r).^2+1/15*(e*r).^3).*exp(-e*r);
% Generate some points in 2 dimensions
%nArr = 1024:1024:10000;
nArr = 12;
m = 50;
nvecs = floor(n/10);
SLQ_reord_K_time = zeros(length(nArr),1);
SLQ_reord_HODLR_time = zeros(length(nArr),1);
SLQ_K_HODLR_time = zeros(length(nArr),1);
SLQ_K_time = zeros(length(nArr),1);

for n=nArr
dim2Points = randn(n,2);
DM = real(DistanceMatrix(dim2Points,dim2Points));
% Create kernel matrix
ep = 10;
K = kernel(ep,DM);
% Bump the eigenvalues
tic
K = K + 1e-1*eye(n);

SLQ_K = SLQ(@(b) K*b,@log,n,m,nvecs);
SLQ_K_time(ii) = toc;
%%
% Preprocessing
% Reordering
tic
reord_K = minDistance2dReordering(dim2Points,kernel,ep);
reord_K = reord_K + 1e-1*eye(n); % keep the reorderings and use as permutation
SLQ_reord_K = SLQ(@(b) reord_K*b,@log,n,m,nvecs);
[sort(eig(K)),sort(eig(reord_K))]
SLQ_reord_K_time(ii) = toc;
approxRank = 10;
diagSize = 130;
I = [1,n];
% maximin reordered kernel matrix (can't do direct kernel b/c eigen-bump.
kMtrxFcn = @(b) reord_K * b;
tic
reord_HODLR = MakeHODLRMtrx(kMtrxFcn,n,approxRank,diagSize,I);
SLQ_reord_HODLR = SLQ(@(b) HODLRMatVec(reord_HODLR,b),@log,n,m,nvecs);
SLQ_reord_HODLR_time(ii) = toc;
% no order kernel matrix
kMtrxFcn = @(b) K * b;
tic
K_HODLR = MakeHODLRMtrx(kMtrxFcn,n,approxRank,diagSize,I);
SLQ_K_HODLR = SLQ(@(b) HODLRMatVec(K_HODLR,b),@log,n,m,nvecs);
SLQ_K_HODLR_time(ii) = toc;
%%
%Accuracy plot
% figure()
% plot(SLQ_K)
% hold on
% plot(SLQ_reord_K)
% plot(SLQ_reord_HODLR)
% plot(SLQ_K_HODLR)
% legend({'SLQ\_K','SLQ\_reord\_K','SLQ\_reord\_HODLR','SLQ\_K\_HOLDR'},'Location','southwest')
% hold off
end
%%
% Speed plot
figure()
plot(SLQ_K_time)
hold on
plot(SLQ_reord_K_time)
plot(SLQ_reord_HODLR_time)
plot(SLQ_K_HODLR_time)
legend({'SLQ\_K','SLQ\_reord\_K','SLQ\_reord\_HODLR','SLQ\_K\_HOLDR'},'Location','southeastoutside')
hold off