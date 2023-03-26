% Matlab code LanczosFullReorthog.m
% For "Applied Numerical Linear Algebra",  Chapter 7, section 3
% Written by James Demmel, Jun  6, 1997
% Extended by Ibrohim Nosirov, Sept 6, 2022
% Perform HODLR-Lanczos with complete reorthogonalization
% the code can be one function.
function T = HODLR_Lanczos(A,q,m)
      qq = q/norm(q);
      Q = qq;
      alpha = zeros(m,1);
      beta = alpha;

      for ii=1:m
            z = HODLRMatVec(A,qq);
            alpha(ii) = qq'*z;
            % reorthogonalize twice to guarantee orthogonality of Lanczos
            % vectors
            z = z - Q*(Q'*z);
            z = z - Q*(Q'*z);
            beta(ii) = norm(z);
            qq = z/beta(ii);
            Q = [Q,qq];
      end
      T = diag(alpha(1:m)) + diag(beta(1:m-1),1) + diag(beta(1:m-1),-1);
end