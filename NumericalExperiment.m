% load a matrix into Disk
% write a read/write anonymous function.
 
close all;
%profile on
nArr = (1:10)*1000;
% make C6 Matern kernel
kernel = @(e,r) (1+e*r+2/5*(e*r).^2+1/15*(e*r).^3).*exp(-e*r); ep=1; d=1;

for ii=1:length(nArr)
    n=nArr(ii); diagSize=130*ii; r=4; I=[1 n^d]; blockSize=500;
    % Point evals at which to sample
    x = CreatePoints(n^d,d,'u');
%    DM = DistanceMatrix(x,x);
%    A = kernel(ep,DM);
%    SaveSquareMatrix(A,blockSize);
    % make kernel matrix and multiply on demand
    kMtrxFcn = @(b) kMtrxMatMat(b,x,kernel,ep,blockSize);
%    base = '/home/ibrohim/Research/S22/peeling_algorithm/';
%    name = ['BigAMatrix',num2str(n),'.h5']
%    fileName = fullfile(base,name);
%    kMtrxFcn = @(b) kMtrxH5MatMat(fileName,b,blockSize,n);
    % sample matrix
    tic
    A = MakeHODLRMtrx(kMtrxFcn,n,r,diagSize,I);
    HODLR_Gamma = HODLR_SLQ(A,n,@log,100,10)
    toc
    tic
    Vanilla_Gamma = SLQ(kMtrxFcn,n,@log,100,10)
    toc
%     DM = DistanceMatrix(x,x);
%     A = kernel(ep,DM);
end

% A = rbf(ep,DM);
% SaveSquareMatrix(A,maxBlockSize); % <- check maxBlockSize
% ReadSquareMatrix('BigAMatrix.h5')
% p = profile("info");
% profsave(p,"myresults3")
% profile off

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

function y = kMtrxH5MatMat(fileName,b,blockSize,n)
    numberOfBlocks = ceil(n/blockSize);
    y = zeros(size(b));
    for j = 1:numberOfBlocks
        startInd = (j-1)*blockSize + 1;
        endInd = min(n,j*blockSize);
        sizeBlock = endInd - startInd+1;

        A_part = h5read(fileName,'/A',[1,startInd],[n,sizeBlock]);
        y = y + A_part*b(startInd:endInd,:);
    end
% read matrix from file to compute
end
