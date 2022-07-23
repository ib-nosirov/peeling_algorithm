% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
n = 1000; diagSize = 125; k = 10; I = [1 n^d];
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
kMtrx = rbf(ep,DM);
output = MakeHODLRMtrx(kMtrx,n^d,k,diagSize,I);
uTree = output(1); zTree= output(2); idxTree = output(3);
% reconstruct the K matrix approximation from U and Z values
kApprox = kMtrx;
kApprox(1:500,501:end) = uTree.get(2) * zTree.get(3)';
kApprox(501:end,1:500) = uTree.get(3) * zTree.get(2)';

kApprox(1:250,251:500) = uTree.get(4) * zTree.get(5)';
kApprox(251:500,1:250) = uTree.get(5) * zTree.get(4)';
kApprox(501:750,751:end) = uTree.get(10) * zTree.get(11)';
kApprox(751:end,501:750) = uTree.get(11) * zTree.get(10)';

kApprox(1:125,126:250) = uTree.get(6) * zTree.get(7)';
kApprox(126:250,1:125) = uTree.get(7) * zTree.get(6)';
kApprox(251:375,376:500) = uTree.get(8) * zTree.get(9)';
kApprox(376:500,251:375) = uTree.get(9) * zTree.get(8)';
kApprox(501:625,626:750) = uTree.get(12) * zTree.get(13)';
kApprox(626:750,501:625) = uTree.get(13) * zTree.get(12)';
kApprox(751:875,876:end) = uTree.get(14) * zTree.get(15)';
kApprox(876:end,751:875) = uTree.get(15) * zTree.get(14)';

norm(kMtrx(1:500,501:end) - uTree.get(2) * zTree.get(3)','fro')...
/norm(kMtrx(1:500,501:end),'fro')
norm(kMtrx(501:end,1:500) - uTree.get(3) * zTree.get(2)','fro')...
/norm(kMtrx(501:end,1:500),'fro')

norm(kMtrx(1:250,251:500) - uTree.get(4) * zTree.get(5)','fro')...
/norm(kMtrx(1:250,251:500),'fro')
norm(kMtrx(251:500,1:250) - uTree.get(5) * zTree.get(4)','fro')...
/norm(kMtrx(251:500,1:250),'fro')
norm(kMtrx(501:750,751:end) - uTree.get(10) * zTree.get(11)','fro')...
/norm(kMtrx(501:750,751:end),'fro')
norm(kMtrx(751:end,501:750) - uTree.get(11) * zTree.get(10)','fro')...
/norm(kMtrx(751:end,501:750),'fro')

norm(kMtrx(1:125,126:250) - uTree.get(6) * zTree.get(7)','fro')...
/norm(kMtrx(1:125,126:250),'fro')
norm(kMtrx(126:250,1:125) - uTree.get(7) * zTree.get(6)','fro')...
/norm(kMtrx(126:250,1:125),'fro')
norm(kMtrx(251:375,376:500) - uTree.get(8) * zTree.get(9)','fro')...
/norm(kMtrx(251:375,376:500),'fro')
norm(kMtrx(376:500,251:375) - uTree.get(9) * zTree.get(8)','fro')...
/norm(kMtrx(376:500,251:375),'fro')
norm(kMtrx(501:625,626:750) - uTree.get(12) * zTree.get(13)','fro')...
/norm(kMtrx(501:625,626:750),'fro')
norm(kMtrx(626:750,501:625) - uTree.get(13) * zTree.get(12)','fro')...
/norm(kMtrx(626:750,501:625),'fro')
norm(kMtrx(751:875,876:end) - uTree.get(14) * zTree.get(15)','fro')...
/norm(kMtrx(751:875,876:end),'fro')
norm(kMtrx(876:end,751:875) - uTree.get(15) * zTree.get(14)','fro')...
/norm(kMtrx(876:end,751:875),'fro')

disp(abs(norm(kMtrx - kApprox, 'fro') / norm(kMtrx, 'fro')))
