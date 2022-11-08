% Matlab code LanczosFullReorthog.m
% For "Applied Numerical Linear Algebra",  Chapter 7, section 3
% Written by James Demmel, Jun  6, 1997
%
% Perform Lanczos with complete reorthogonalization
function T = HODLR_Lanczos(A,q,m)
      qq = q/norm(q);
      Q = qq;
      alpha = zeros(m,1);
      beta = alpha;

      for i=1:m
            z = HODLRMatVec(A,qq);
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