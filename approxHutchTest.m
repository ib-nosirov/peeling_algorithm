% Construct Kernel matrix
N = 128;
K = zeros(N);
sigma = N/10;
for ii = 1:N
    for jj = 1:N
        K(ii,jj) = exp(-(ii-jj)^2/sigma^2);
    end
end
figure(1);clf;imagesc(K);colorbar;title('Kernel matrix K');

% Use Martinson's technique to approximate the top right and bottom left
% quadrants of the Kernel matrix. Actually we just use the "U" part; we
% don't seem to need the Z part here for a one-level decomposition.
r = 5;
Omega = zeros(N,2*r);
Omega(1:N/2,r+1:2*r) = randn(N/2,r);
Omega(N/2+1:N,1:r) = randn(N/2,r);
Y = K*Omega; % THIS USES 2r MULTIPLIES BY K

Y2 = Y(1:N/2,1:r); 
U2 = orth(Y2);
Y3 = Y(N/2+1:N,r+1:2*r);
U3 = orth(Y3);

% Stack the orthobases
Q = [U2 zeros(N/2,r); zeros(N/2,r) U3];
figure(2);clf;imagesc(Q);colorbar;title('Approx. basis matrix Q');

% Confirm that tr(K) = tr(Q'*K*Q) + (I-Q*Q')*K*(I-Q*Q')
abs(trace(K) - (trace(Q'*K*Q) + trace((eye(N)-Q*Q')*K*(eye(N)-Q*Q'))))

% Part 1: Compute tr(Q'*K*Q) exactly
tracePart1 = trace(Q'*K*Q) % THIS USES 2r MULTIPLIES BY K

% Part 2: Estimate tr((I-Q*Q')*K*(I-Q*Q')) via Hutchinson's Estimator
G = randn(N,2*r); % can change the number of columns here
tracePart2 = (1/size(G,2))*trace(G'*(eye(N)-Q*Q')*K*(eye(N)-Q*Q')*G) % THIS USES 2r MULTIPLIES BY K

trace(K)
traceEst = tracePart1 + tracePart2

% The variance of this estimator should be proportional to
% ||((I-Q*Q')*K*(I-Q*Q'))||_F^2 which appears to be much smaller than
% tr(K)^2
varOriginal = (trace(K))^2
varEst = norm((eye(N)-Q*Q')*K*(eye(N)-Q*Q'),'fro')^2
figure(3);clf;imagesc((eye(N)-Q*Q')*K*(eye(N)-Q*Q'));colorbar;title('Residual (I-Q*Q'')*K*(I-Q*Q'')');

figure(4);clf;imagesc(Q*Q'*K*Q*Q');title('Q*Q''*K*Q*Q''');

figure(5);clf;plot(svd(K));

