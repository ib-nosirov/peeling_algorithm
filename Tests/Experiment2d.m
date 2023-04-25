% Gaussian kernel
% mirror lower triangular to upper triangular (is this necessary?)
% Test MakeHODLRMtrx
close all;
nArr = [1024,2048,4096];
for ii=1:length(nArr)
	n = nArr(ii); d=1; diagSize=130; r=10; I=[1 n];
    dim2Points = randn(n,2);
    radiusVec = dim2Points - dim2Points;
    rbfKernel = @(x,y) exp(-norm(x-y)^2);
    nonReorderMtrx = zeros(n);
    for kk=1:n
        for jj=1:n
            nonReorderMtrx(kk,jj) = rbfKernel(dim2Points(kk,:),dim2Points(jj,:));
        end
    end
    semilogy(svd(nonReorderMtrx))
    imagesc(nonReorderMtrx)
    %%
    M = minDistance2dReordering(dim2Points,rbfKernel);
    M = M./(1e-6 + vecnorm(M));
    %M = M + rot90(M,2);
    %M = tril(M)+tril(M,-1)';
    %M = M - mean(M(:));
    imagesc(M)

    %$imagesc(M);
    kMtrxFcn = @(b) M*b;
	% Point evals at which to sample
	%x = CreatePoints(n^d,d,'u');
	% compute the absolute difference
	%DM = DistanceMatrix(x,x);
	% sample matrix
	%kMtrx = kMtrxFcn(ep,DM);

    % make a cell array of the diagonals of M.


	K = MakeHODLRMtrx(kMtrxFcn,n,r,diagSize,I);

	% reconstruct the K matrix approximation from U and Z values
	kApprox = zeros(n,n);
	uTree = K{1};
	zTree = K{2};
	leavesCell = K{3};
	idxTree = K{4};
	it = idxTree.breadthfirstiterator;
	treeDepth = floor(log2(nnodes(idxTree)+1));
	offset = length(it)-2^(treeDepth-1);
figure()
    for idx=3:2:nnodes(idxTree)-2^(treeDepth-1)
		s1 = table(idxTree.get(it(idx-1))).Var1(1);
		f1 = table(idxTree.get(it(idx-1))).Var1(2);
		s2 = table(idxTree.get(it(idx))).Var1(1);
		f2 = table(idxTree.get(it(idx))).Var1(2);
		kApprox(s1:f1,s2:f2) = uTree.get(it(idx-1)) * zTree.get(it(idx))';
		kApprox(s2:f2,s1:f1) = uTree.get(it(idx)) * zTree.get(it(idx-1))';
        %% testing
        tmp = M(s1:f1,s2:f2);
        tmp = tmp - mean(tmp(:));
        semilogy(svd(tmp))
        hold on
        %plot(real(fft(tmp)))
        % 2d fft then vectorize, or compute singular values.
        tmp = M(s2:f2,s1:f1);
        tmp = tmp - mean(tmp(:));
        semilogy(svd(tmp))
keyboard
        %plot(real(fft(tmp(:))))
    end
    for idx=1:length(leavesCell)
		s = table(idxTree.get(it(idx+offset))).Var1(1);
		f = table(idxTree.get(it(idx+offset))).Var1(2);
		kApprox(s:f,s:f) = leavesCell{idx};
        %figure()
        tmp = M(s:f,s:f);
        %semilogy(svd(tmp))
        %plot(real(fft(tmp(:))))
    end
    %imagesc(kApprox-M)
    %colorbar
    keyboard
    % reconstruct matrix by matrix multiply
    kApproxMatVec = HODLRMatVec(K,eye(n));
    imagesc(kApproxMatVec)
	imagesc(M-kApproxMatVec)

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
	
    % Test HODLRMatVec
    b = randn(n,1);
    yApprox = HODLRMatVec(K,b);
    yExact = M*b;
    matVecRelativeError = abs(norm(yExact-yApprox,'fro')/norm(yExact,'fro'));

    % Test Lanczos and HODLR-Lanczos % Relative error of 0.4007
    q = randn(n,1);
    m = 10;
    exactLanczos = Lanczos(kMtrxFcn,q,m);
    approxLanczos = Lanczos(@(b) HODLRMatVec(K,b),q,m);
    relativeLanczosError = abs(norm(exactLanczos-approxLanczos,'fro')/norm(exactLanczos,'fro'));
    
    % Test HODLR-SLQ and SLQ
    MATLAB_Gamma = trace(logm(M));
    HODLR_Gamma = SLQ(@(b) HODLRMatVec(K,b),@log,n,50,100);
    SLQ_Gamma = SLQ(kMtrxFcn,@log,n,50,100);
end

% compare against MATLAB's reorderings.
% run the numerical experiments
% for 1d and 2d, bump the eigenvalues.
