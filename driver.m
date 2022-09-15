close all;
nArr = [400,800,1600,3200];
%for ii=1:length(nArr)
% make a kernel
rbf = @(e,r) exp(-e*r.^2); ep=100; d=1;
n = 1000; diagSize=125; k=10; I=[1 n^d]; maxBlockSize = 100;
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
kMtrxFcn = @(b) kMtrxMatVec(b,x,rbf,ep,maxBlockSize);

%A = MakeHODLRMtrx(kMtrxFcn,n^d,k,diagSize,I);
b = ones(n,1);
yApprox = kMtrxMatVec(b,x,rbf,ep,maxBlockSize);
%yApprox = HODLRMatVec(A,b);
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
kMtrx = rbf(ep,DM);
yExact = kMtrx * b;
norm(yExact-yApprox)

function y = kMtrxMatVec(b,x,rbf,ep,maxBlockSize)
    y = zeros(size(b));
    indices = 1:maxBlockSize:length(y)+1;
    for ii = 2:length(indices)
        iIdx = indices(ii)-1;
        iPrevIdx = indices(ii-1);
        for jj = 2:length(indices)
            jIdx = indices(jj)-1;
            jPrevIdx = indices(jj-1);
            kMtrxBlock = rbf(ep,DistanceMatrix(x(iPrevIdx:iIdx),x(jPrevIdx:jIdx)));
            y(iPrevIdx:iIdx,:) = y(iPrevIdx:iIdx,:) + kMtrxBlock*b(iPrevIdx:iIdx,:);
        end
    end
end
% reconstruct the K matrix approximation from U and Z values
% kApproxDims = size(kMtrx);
% kApprox = zeros(kApproxDims(1),kApproxDims(2));
% it = idxTree.breadthfirstiterator;
% treeDepth = floor(log2(nnodes(idxTree)+1));
% offset = length(it)-2^(treeDepth-1);
% 
% for idx=3:2:nnodes(idxTree)-2^(treeDepth-1)
% 	s1 = table(idxTree.get(it(idx-1))).Var1(1);
% 	f1 = table(idxTree.get(it(idx-1))).Var1(2);
% 	s2 = table(idxTree.get(it(idx))).Var1(1);
% 	f2 = table(idxTree.get(it(idx))).Var1(2);
% 	kApprox(s1:f1,s2:f2) = uTree.get(it(idx-1)) * zTree.get(it(idx))';
% 	kApprox(s2:f2,s1:f1) = uTree.get(it(idx)) * zTree.get(it(idx-1))';
% end
% for idx=1:length(leavesCell)
% 	s = table(idxTree.get(it(idx+offset))).Var1(1);
% 	f = table(idxTree.get(it(idx+offset))).Var1(2);
% 	kApprox(s:f,s:f) = leavesCell{idx};
% end
% disp(abs(norm(kMtrx-kApprox,'fro') / norm(kMtrx,'fro')))
% figure(ii)
% imagesc(kMtrx-kApprox)
% end