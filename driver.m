% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
n = 1000; diagSize = 125; k = 10; I = [1 n^d];
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
kMtrx = rbf(ep,DM);
output = MakeHODLRMtrx(kMtrx,n^d,k,diagSize,I);
uTree = output{1};
zTree= output{2};
leavesCell = output{3};
idxTree = output{4};
% reconstruct the K matrix approximation from U and Z values
kApproxDims = size(kMtrx);
kApprox = zeros(kApproxDims(1),kApproxDims(2));
it = idxTree.breadthfirstiterator;
treeDepth = floor(log2(nnodes(idxTree)+1));
offset = length(it)-2^(treeDepth-1)

for idx=3:2:nnodes(idxTree)-2^(treeDepth-1)
	s1 = table(idxTree.get(it(idx-1))).Var1(1);
	f1 = table(idxTree.get(it(idx-1))).Var1(2);
	s2 = table(idxTree.get(it(idx))).Var1(1);
	f2 = table(idxTree.get(it(idx))).Var1(2);
	kApprox(s1:f1,s2:f2) = uTree.get(it(idx-1)) * zTree.get(it(idx))';
	kApprox(s2:f2,s1:f1) = uTree.get(it(idx)) * zTree.get(it(idx-1))';
end
for idx=1:length(leavesCell)
	s = table(idxTree.get(it(idx+offset))).Var1(1);
	f = table(idxTree.get(it(idx+offset))).Var1(2);
	kApprox(s:f,s:f) = leavesCell{idx};
end
disp(abs(norm(kMtrx - kApprox, 'fro') / norm(kMtrx, 'fro')))
imagesc(kMtrx - kApprox)
