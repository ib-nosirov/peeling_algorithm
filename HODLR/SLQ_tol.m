function [gamma,ii] = SLQ_tol(linearOperator,f,n,m,tol,groundTruth)
%% function ld = Lanc_Quad_LogDet( A, m, nvecs)
% The function computes an approximate log-determinant of PSD matrix A
% using the Stochastic Lanczos Quadrature Approximation
%-- Inputs
% linearOperator - a hodlr matrix multiply function
% n = 
% m - Number of Lanczos steps (degree)
% nvecs - Number of starting vectors
%-- Output
% ld - The logdeterminant of A estimated by SLQ
% z1 - Individual estimates for each starting vector v_l


	%% Initialization
	cntEst = 0;
	%% Main loop
    ii = 1;
    gamma = 0;
	while abs(groundTruth - gamma) > tol && ii < 3*n/10
        if mod(ii,n/10) == 0 
            cnEst = 0;
            ii = 1;
            gamma = 0;
        end
		w = sign(randn(n,1)); % Random rademacher vector
		v0 = w /norm(w);
		[H,V,g] = Lanczos(linearOperator,v0,m); % m steps of Lanczos algorithm
		%H = (H+H')/2;
		[eigvec,D]=eig(H);
		theta  = abs(diag(D));
		gamma = eigvec(1,:).^2;
		% sum of gamma2*log(theta)
		count = sum(gamma'.*(f(theta)));
		cntEst = (cntEst+count);
		gamma = n*(cntEst/ii);
		ii = ii+1;
    end
end