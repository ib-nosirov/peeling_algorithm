
% choose a set of matrix sizes (make sure to include prime numbers)
nArr = 10:12; % 2^16 is too large to store on RAM.
nArr = 2.^nArr;
SLQ_data = nArr;
HODLRSLQ_data = nArr;
MATLAB_data = nArr;
jj = 1;
m = 10;
d = 130;
r = 10;
ep=10;
for n=nArr
    % fix a compression rate by fixing the size of the main diagonal blocks
    % approximate rank should, in theory, should stay contant or grow
    % logarithmically.

    I=[1,n];
    % Point evals at which to sample
    x = CreatePoints(n,1,'u');
    % compute the absolute difference
    DM = DistanceMatrix(x,x);
    % C6 Matern.
    kernel = @(e,r) (1+e*r+2/5*(e*r).^2+1/15*(e*r).^3).*exp(-e*r);
    M = kernel(ep,DM);
    M = M + 1e-1*eye(n);
    %kMtrxFcn = @(b) kernel(ep,DistanceMatrix(x,b)) + 1e-1*b; % issue: b is not
    %always a vector; sometimes a matrix.
    kMtrxFcn = @(b) M*b;

    % reconstruct this.
    %nvecs = floor(n/10);
    groundTruth = trace(logm(M));
    tic
    [SLQ_Gamma,numIters] = SLQ_tol(@(b)M*b,@log,n,m,1e-1,groundTruth);
    SLQ_data(jj) = toc;
    numIters
    SLQ_error = abs(trace(logm(M))-SLQ_Gamma)
    
    tic
    K = MakeHODLRMtrx(kMtrxFcn,n,r,d,I);
    [HODLR_Gamma,numIters] = SLQ_tol(@(b) HODLRMatVec(K,b),@log,n,m,1e-1,groundTruth);
    HODLRSLQ_data(jj) = toc;
    numIters
    HODLRSLQ_error = abs(trace(logm(M))-HODLR_Gamma)
    
    tic
    matlabGamma = trace(logm(M));
    MATLAB_data(jj) = toc;
    
    jj = jj+1;
end
plot(SLQ_data)
hold on
plot(HODLRSLQ_data)
plot(MATLAB_data)
legend('SLQ\_data','HODLR\_SLQ\_data','MATLAB')
hold off