close all;
profile on
nArr = [400,800,1600,3200];
%for ii=1:length(nArr)
% make a kernel
rbf = @(e,r) exp(-e*r.^2); ep=6000; d=1;
n = 6400; diagSize=125; k=10; I=[1 n^d]; maxBlockSize = 800;
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
kMtrxFcn = @(b) kMtrxMatMat(b,x,rbf,ep,maxBlockSize);

% sample matrix
A = MakeHODLRMtrx(kMtrxFcn,n,k,diagSize,I);
Gamma = HODLR_SLQ(A,n,@log,20,10);
DM = DistanceMatrix(x,x);
A = rbf(ep,DM);
%A = orth(randn(n)); A = A*diag(n+1:2*n)*A';
SLQ(A,n,@log,20,10);
p = profile("info");
profsave(p,"myresults3")
profile off
function y = kMtrxMatMat(b,x,rbf,ep,maxBlockSize)
    y = zeros(size(b));
    indices = 1:maxBlockSize:length(y)+1;
    for ii = 2:length(indices)
        iIdx = indices(ii)-1;
        iPrevIdx = indices(ii-1);
        kMtrxBlock = rbf(ep,DistanceMatrix(x,x(iPrevIdx:iIdx)));
        y = y + kMtrxBlock*b(iPrevIdx:iIdx,:);
    end
end
