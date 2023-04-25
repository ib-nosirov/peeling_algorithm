function [ld,z1] = SLQ(linearOperator,f,n,m,nvecs)
%% function ld = Lanc_Quad_LogDet( A, m, nvecs)
% The function computes an approximate log-determinant of PSD matrix A
% using the Stochastic Lanczos Quadrature Approximation
%-- Inputs
% A - the symmetric positive definite input matrix
% m - Number of Lanczos steps (degree)
% nvecs - Number of starting vectors
%-- Output
% ld - The logdeterminant of A estimated by SLQ
% z1 - Individual estimates for each starting vector v_l


	%% Initialization
	cnt_est=0;
	%% Main loop
	for ii = 1:nvecs
        w = double(1:n == ii)';
		%w = sign(randn(n,1)); % Random radamacher vector
		v0 = w /norm(w);
		[H,V,g] = Lanczos(linearOperator,v0,m); % m steps of Lanczos algorithm
		%H = (H+H')/2;
		[eigvec,D]=eig(H);
		theta  = abs(diag(D));
		gamma2 = eigvec(1,:).^2;
		%% sum of gamma2*log(theta)
		count=sum(gamma2'.*(f(theta)));
		z1(ii) = (count)*n;
		cnt_est = (cnt_est+count);
		zz(ii) = n*(cnt_est/ii);
	end
	ld=zz;
