% Matlab code LanczosFullReorthog.m
% For "Applied Numerical Linear Algebra",  Chapter 7, section 3
% Written by James Demmel, Jun  6, 1997
%
% Perform Lanczos with complete reorthogonalization
%
% Inputs:
%
%   A = input matrix (may be sparse)
%   e = true eigenvalues, sorted from largest to smallest
%   q = starting vector for Lanczos
%   m = number of Lanczos steps to take
%   mm1,mm2 = indices of extreme eigenvalues for which to plot 
%             error and error bounds
%   mine, maxe = range in which to plot eigenvalues
%   steplabel = array of values at which to draw vertical lines in plots
%   recomp = 1 to perform algorithm from scratch
%          = 0 just to redisplay data, with possibly different
%               mine, maxe, mm1, mm2, steplabel
%   whichbound = 2 to plot all errors and error bound
%              = 1 to plot local error and error bounds
%              = 0 to plot error bounds only
%
% Outputs:
%
%   plot of Ritz values for all m Lanczos steps, in range [mine,maxe]
%
%   plot of global error (distance from final eigenvalue) (if whichbound=2)
%           local error (distance from closest eigenvalue) (if whichbound=1)
%           error bound (provided by Lanczos algorithm)
%     for mm1 to mm2 smallest and mm1 to mm2 largest eigenvalues
%
%   plot of components of each eigenvector in current Lanczos vector
%           (only valid if A is diagonal with eigenvalues sorted
%            from largest to smallest!)
%     for mm1 to mm2 smallest and mm1 to mm2 largest eigenvalues
% 
% Lanczos needs to take in function handle rather than function

function T = Lanczos(kMtrxFcn,q,m)
      qq = q/norm(q);
      Q = qq;
      alpha = zeros(m,1);
      beta = alpha;

      for i=1:m
            if (rem(i,10)==0)
            end
            z = kMtrxFcn(qq);
            alpha(i) = qq'*z;
            % reorthogonalize twice to guarantee orthogonality of Lanczos
            % vectors
            z = z - Q*(Q'*z);
            z = z - Q*(Q'*z);
            beta(i) = norm(z);
            qq = z/beta(i);
            Q = [Q,qq];
      end
      T = diag(alpha(1:m)) + diag(beta(1:m-1),1) + diag(beta(1:m-1),-1);
end
