
% choose a set of matrix sizes (make sure to include prime numbers)
nArr = 10:15; % 2^16 is too large to store on RAM.
nArr = 2.^nArr;
for n=nArr
    % fix a compression rate by fixing the size of the main diagonal blocks
    d = floor(n/10);
    % approximate rank should, in theory, also grow with matrix size.
    r = floor(d/10);
    ep=10;
    I=[1,n];
    % Point evals at which to sample
    x = CreatePoints(n,1,'u');
    % compute the absolute difference
    DM = DistanceMatrix(x,x);
    % C6 Matern.
    kernel = @(e,r) (1+e*r+2/5*(e*r).^2+1/15*(e*r).^3).*exp(-e*r);
    M = kernel(ep,DM);
    M = M + 1e-1*eye(n);
    kMtrxFcn = @(b) kernel(ep,x-b) + 1e-1*b; % kernel(something)
    % reconstruct this.
    nvecs = floor(n/10);
    m = floor(nvecs/5);
    tic
    SLQ_Gamma = SLQ(@(b)M*b,@log,n,m,nvecs);
    toc
    tic
    K = MakeHODLRMtrx(kMtrxFcn,n,r,d,I);
    HODLR_Gamma = SLQ(@(b) HODLRMatVec(K,b),@log,n,m,nvecs);
    toc
    tic
    matlabGamma = trace(logm(M));
    toc
    trace(logm(M))-HODLR_Gamma(end)
end