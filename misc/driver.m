% load a matrix into Disk
% write a read/write anonymous function.
 
close all;
%profile on
%nArr = (1:10)*1000;
% make C6 Matern kernel
%kernel = @(e,r) (1+e*r+2/5*(e*r).^2+1/15*(e*r).^3).*exp(-e*r); ep=10; d=1;
%for ii=1:length(nArr)
d = 1;
    n=1000; diagSize=130*1; r=4; I=[1,n^d]; maxBlockSize=500;
    M = randn(n); M = M*M';
    % Point evals at which to sample
    x = CreatePoints(n^d,d,'u');
    % make kernel matrix and multiply on demand
    kMtrxFcn = @(b) M*b;
    % sample matrix; has to be logm here since the below 2 lines are
    % equivalent.
    trace(logm(M))
    %sum(log(eig(M)))
%    tic
    M = MakeHODLRMtrx(kMtrxFcn,n,r,diagSize,I);
    % has to be @log here since we take log(eig(T))
    HODLR_Gamma = HODLR_SLQ(M,n,@log,50,10)
    Gamma = SLQ(kMtrxFcn,n,@log,50,10)
%    toc
%end

function y = kMtrxMatMat(b,x,kernel,ep,maxBlockSize)
    y = zeros(size(b));
    indices = 1:maxBlockSize:length(y)+1;
    for ii = 2:length(indices)
        iIdx = indices(ii)-1;
        iPrevIdx = indices(ii-1);
        kMtrxBlock = kernel(ep,DistanceMatrix(x,x(iPrevIdx:iIdx)));
        y = y + kMtrxBlock*b(iPrevIdx:iIdx,:);
    end
end