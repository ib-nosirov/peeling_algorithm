load('thermomech_TC.mat')
A = Problem.A;
%A = A*A';
kMtrxFcn = @(b) A*b;
n=102158;
d=1;
k=5;
diagSize=150;
m=50;
nvecs=100;
I=[1 n];
%%
K = MakeHODLRMtrx(kMtrxFcn,n,k,diagSize,I);
% 	% reconstruct the K matrix approximation from U and Z values
% 	kApprox = zeros(n,n);
% 	uTree = K{1};
% 	zTree = K{2};
% 	leavesCell = K{3};
% 	idxTree = K{4};
% 	it = idxTree.breadthfirstiterator;
% 	treeDepth = floor(log2(nnodes(idxTree)+1));
% 	offset = length(it)-2^(treeDepth-1);
%     % this code works. relative reconstruction error = 4.4765e-14
%     for idx=3:2:nnodes(idxTree)-2^(treeDepth-1)
% 		s1 = table(idxTree.get(it(idx-1))).Var1(1);
% 		f1 = table(idxTree.get(it(idx-1))).Var1(2);
% 		s2 = table(idxTree.get(it(idx))).Var1(1);
% 		f2 = table(idxTree.get(it(idx))).Var1(2);
% 		kApprox(s1:f1,s2:f2) = uTree.get(it(idx-1)) * zTree.get(it(idx))';
% 		kApprox(s2:f2,s1:f1) = uTree.get(it(idx)) * zTree.get(it(idx-1))';
%     end
%     for idx=1:length(leavesCell)
% 		s = table(idxTree.get(it(idx+offset))).Var1(1);
% 		f = table(idxTree.get(it(idx+offset))).Var1(2);
% 		kApprox(s:f,s:f) = leavesCell{idx};
%     end
%     % reconstruct matrix by matrix multiply
%     kApproxMatVec = HODLRMatVec(K,eye(n));
%     norm(A-kApproxMatVec,'fro')/norm(A,'fro')
% %MATLAB_Gamma = trace(logm(A))
%%
% rng(1);
%[ldUbaru,z1] = Lanc_Quad_LogDet(A,m,nvecs);
[hodlr_ld,z1] = SLQ(@(b) HODLRMatVec(K,b),n,m,nvecs);
%[ld,z1] = SLQ(kMtrxFcn,@log,n,m,nvecs);
figure()
hold on
for ii=1:3
[ldUbaru,z1] = Lanc_Quad_LogDet(A,m,nvecs);
plot(ldUbaru)
end
%plot(hodlr_ld)
%plot(ld)
%plot(ldUbaru)
%legend({'SLQ','Ubaru'},'Location','northeast')
hold off
ldUbaru(end)
ld(end)
%hodlr_ld(end)
