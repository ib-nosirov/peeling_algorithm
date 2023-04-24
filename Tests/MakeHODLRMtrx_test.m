function MakeHODLR_test()
n = nArr(ii); d=1; diagSize=102158; r=5; I=[1 n^d]; % reconstruction error is governed by rank.
	%kMtrxFcn = @(b) M*b;
    
	% Point evals at which to sample
	x = CreatePoints(n^d,d,'u');
	% compute the absolute difference
	DM = DistanceMatrix(x,x);
    kernel = @(e,r) (1+e*r+2/5*(e*r).^2+1/15*(e*r).^3).*exp(-e*r); ep=10; d=1;
    M = kernel(ep,DM);
    M = M + 1e-1*eye(n);
%    M = M + triu(M);
    kMtrxFcn = @(b) M*b;
	% sample matrix
	%kMtrx = kMtrxFcn(ep,DM);
    % Assume it's symmetric
	K = MakeHODLRMtrx(kMtrxFcn,n^d,r,diagSize,I);

	% reconstruct the K matrix approximation from U and Z values
	kApprox = zeros(n,n);
	uTree = K{1};
	zTree = K{2};
	leavesCell = K{3};
	idxTree = K{4};
	it = idxTree.breadthfirstiterator;
	treeDepth = floor(log2(nnodes(idxTree)+1));
	offset = length(it)-2^(treeDepth-1);
    % this code works. relative reconstruction error = 4.4765e-14
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
    % reconstruct matrix by matrix multiply
    kApproxMatVec = HODLRMatVec(K,eye(n));
end