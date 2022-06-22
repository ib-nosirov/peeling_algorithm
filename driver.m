n = 1000; diag_size = 125; k = 10; I = [1 n];
% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
K_mtrx = rbf(ep,DM);
[U_tree Z_tree idx_tree f_approx] = NumericalExperiment(K_mtrx,n,k,diag_size,I);

% reconstruct the K matrix approximation from U and Z values
%K_approx = K_mtrx;
%K_approx(1:500,501:end) = U_tree.get(2) * Z_tree.get(3)';
%K_approx(501:end,1:500) = U_tree.get(3) * Z_tree.get(2)';
%
%K_approx(1:250,251:500) = U_tree.get(4) * Z_tree.get(5)';
%K_approx(251:500,1:250) = U_tree.get(5) * Z_tree.get(4)';
%K_approx(501:750,751:end) = U_tree.get(10) * Z_tree.get(11)';
%K_approx(751:end,501:750) = U_tree.get(11) * Z_tree.get(10)';
%
%K_approx(1:125,126:250) = U_tree.get(6) * Z_tree.get(7)';
%K_approx(126:250,1:125) = U_tree.get(7) * Z_tree.get(6)';
%K_approx(251:375,376:500) = U_tree.get(8) * Z_tree.get(9)';
%K_approx(376:500,251:375) = U_tree.get(9) * Z_tree.get(8)';
%K_approx(501:625,626:750) = U_tree.get(12) * Z_tree.get(13)';
%K_approx(626:750,501:625) = U_tree.get(13) * Z_tree.get(12)';
%K_approx(751:875,876:end) = U_tree.get(14) * Z_tree.get(15)';
%K_approx(876:end,751:875) = U_tree.get(15) * Z_tree.get(14)';

%norm(K_mtrx(1:500,501:end) - U_tree.get(2) * Z_tree.get(3)','fro')...
%/norm(K_mtrx(1:500,501:end),'fro')
%norm(K_mtrx(501:end,1:500) - U_tree.get(3) * Z_tree.get(2)','fro')...
%/norm(K_mtrx(501:end,1:500),'fro')
%
norm(K_mtrx(1:250,251:500) - U_tree.get(4) * Z_tree.get(5)','fro')...
/norm(K_mtrx(1:250,251:500),'fro')
norm(K_mtrx(251:500,1:250) - U_tree.get(5) * Z_tree.get(4)','fro')...
/norm(K_mtrx(251:500,1:250),'fro')
%norm(K_mtrx(501:750,751:end) - U_tree.get(10) * Z_tree.get(11)','fro')...
%/norm(K_mtrx(501:750,751:end),'fro')
%norm(K_mtrx(751:end,501:750) - U_tree.get(11) * Z_tree.get(10)','fro')...
%/norm(K_mtrx(751:end,501:750),'fro')
%
%norm(K_mtrx(1:125,126:250) - U_tree.get(6) * Z_tree.get(7)','fro')...
%/norm(K_mtrx(1:125,126:250),'fro')
%norm(K_mtrx(126:250,1:125) - U_tree.get(7) * Z_tree.get(6)','fro')...
%/norm(K_mtrx(126:250,1:125),'fro')
%norm(K_mtrx(251:375,376:500) - U_tree.get(8) * Z_tree.get(9)','fro')...
%/norm(K_mtrx(251:375,376:500),'fro')
%norm(K_mtrx(376:500,251:375) - U_tree.get(9) * Z_tree.get(8)','fro')...
%/norm(K_mtrx(376:500,251:375),'fro')
%norm(K_mtrx(501:625,626:750) - U_tree.get(12) * Z_tree.get(13)','fro')...
%/norm(K_mtrx(501:625,626:750),'fro')
%norm(K_mtrx(626:750,501:625) - U_tree.get(13) * Z_tree.get(12)','fro')...
%/norm(K_mtrx(626:750,501:625),'fro')
%norm(K_mtrx(751:875,876:end) - U_tree.get(14) * Z_tree.get(15)','fro')...
%/norm(K_mtrx(751:875,876:end),'fro')
%norm(K_mtrx(876:end,751:875) - U_tree.get(15) * Z_tree.get(14)','fro')...
%/norm(K_mtrx(876:end,751:875),'fro')
%
%imagesc(K_mtrx - K_approx)
%disp(abs(norm(K_mtrx - K_approx, 'fro') / norm(K_mtrx, 'fro')))
%%end
