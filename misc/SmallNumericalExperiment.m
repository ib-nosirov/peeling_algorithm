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
    %x = CreatePoints(n^d,d,'u');
    % make kernel matrix and multiply on demand
    kMtrxFcn = @(b) M*b;
    % sample matrix
%    tic
    A = MakeHODLRMtrx(kMtrxFcn,n,r,diagSize,I);
    HODLR_Gamma = HODLR_SLQ(A,n,@log,50,10)
    Vanilla_Gamma = SLQ(kMtrxFcn,n,@log,100,10)
    trace(log(M))
%    toc
%end

function y = kMtrxMatMat(b,x,kernel,ep,blockSize)
    y = zeros(size(b));
    indices = 1:blockSize:length(y)+1;
    for ii = 2:length(indices)
        iIdx = indices(ii)-1;
        iPrevIdx = indices(ii-1);
        kMtrxBlock = kernel(ep,DistanceMatrix(x,x(iPrevIdx:iIdx)));
        y = y + kMtrxBlock*b(iPrevIdx:iIdx,:);
    end
end