% Gaussian kernel
% mirror lower triangular to upper triangular (is this necessary?)
% Test MakeHODLRMtrx
close all;
nArr = [1024,2048,4096];
for ii=1:length(nArr)
	n = nArr(ii); d=2; diagSize=130; r=10; I=[1 n^d];
    dim2Points = randn(n,2);
    rbfKernel = @(x,y) exp(-norm(x-y)^2);
    M = minDistance2dReordering(dim2Points,rbfKernel);
    M = tril(M)+tril(M,-1)';
    kMtrxFcn = @(b) M*b;
	% Point evals at which to sample
	%x = CreatePoints(n^d,d,'u');
	% compute the absolute difference
	%DM = DistanceMatrix(x,x);
	% sample matrix
	%kMtrx = kMtrxFcn(ep,DM);

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
	relativeReconstructionError = abs(norm(M-kApprox,'fro')/norm(M,'fro'))
	figure(ii)
	imagesc(M-kApprox)
    
    % Test HODLRMatVec
    b = ones(n,1);
    yApprox = HODLRMatVec(K,b);
    yExact = M*b;
    matVecRelativeError = abs(norm(yExact-yApprox,'fro')/norm(yExact,'fro'))

    % Test Lanczos and HODLR-Lanczos % Relative error of 0.4007
    q = randn(n,1);
    m = 10;
    exactLanczos = Lanczos(kMtrxFcn,q,m);
    approxLanczos = HODLR_Lanczos(K,q,m);
    relativeLanczosError = abs(norm(exactLanczos-approxLanczos,'fro')/norm(exactLanczos,'fro'))
    
    % Test HODLR-SLQ and SLQ
    MATLAB_Gamma = trace(logm(M))
    HODLR_Gamma = HODLR_SLQ(K,n,@log,50,100)
    SLQ_Gamma = SLQ(kMtrxFcn,n,@log,50,100)
end