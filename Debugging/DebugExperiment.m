% Test MakeHODLRMtrx
close all;
nArr = [128]; % increase this for speed test.
mArr = 10; % we start seeing 2 orders of magnitude dropoff when m < 300.
nvArr= 100; % not much of a difference, nv > m, but m should also be large.
numTrials = 5;
accuracy = zeros(length(mArr),length(nvArr));
% need more samples to know for sure.
data = zeros(length(nArr),length(mArr),length(nvArr));
for ii=1:length(nArr)
	n = nArr(ii); d=1; diagSize=16; r=4; I=[1 n^d]; % reconstruction error is governed by rank.
    fprintf('Matrix size %d\n',n)
	M = randn(n,10);
    M = M*M';
    M = M + 1e-3*eye(n);
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

	%relativeReconstructionError = abs(norm(M-kApprox,'fro')/norm(M,'fro'))
	%figure(ii)
	%imagesc(abs(M-kApprox)./abs(M))
    %colorbar
    
    % Test HODLRMatVec    
    b = randn(n,1);
    yApprox = HODLRMatVec(K,b);
    yExact = M*b;
    matVecRelativeError = abs(norm(yExact-yApprox,'fro')/ ...
        norm(yExact,'fro'));
%%
    % Test Lanczos and HODLR-Lanczos % Relative error of 1.00e-14.
    q = randn(n,1);
    m = 10;
    exactLanczos = Lanczos(kMtrxFcn,q,m);
    approxLanczos = Lanczos(@(b) HODLRMatVec(K,b),q,m);
    relativeLanczosError = abs(norm(exactLanczos-approxLanczos,'fro') ...
        /norm(exactLanczos,'fro'));
 %%
    MATLAB_Gamma = trace(logm(M));
    for jj=1:length(mArr)
        for kk=1:length(nvArr)
            trialsSum = 0;
            for ll=1:numTrials
            % Test HODLR-SLQ and SLQ
            %tic
            %MATLAB_Gamma = trace(logm(M));
            %toc
            %tic
            %SLQ_Gamma = SLQ(kMtrxFcn,n,@log,mArr(jj),nvArr(kk));
            %toc
            %MATLAB_Gamma-SLQ_Gamma
            %tic
            HODLR_Gamma = SLQ(@(b) HODLRMatVec(K,b),@log,n,mArr(jj), ...
            nvArr(kk));
            %toc
            trialsSum = trialsSum + abs(MATLAB_Gamma-HODLR_Gamma);
            end
            accuracy(jj,kk) = table(trialsSum/numTrials).Var1(end);
        end
    end
end
figure(1)
imagesc(accuracy)
% to check when variability goes down, we will start at 500,500 and go
% down sequentially. We will also change matrix size.

% I need to show, through my numerical experiments that we can't do much
% better than within 10.

    % TODO:
    % 6. fix the zero-padding issues (U matrix inside of MakeHODLRMatrix)
    % 8. Functions of Matrices J. Higham (f(A)b problem)
    % TODO:
    % 1. Try to replicate Ubaru results.