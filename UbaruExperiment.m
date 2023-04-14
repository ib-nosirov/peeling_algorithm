load('../../data/California.mat')
A = Problem.A;
kMtrxFcn = @(b) A*b;
n=9664;
d=1;
r=5;
diagSize=151;
m=55;
nvecs=100;
I=[1 n^d];
K = MakeHODLRMtrx(kMtrxFcn,n^d,r,diagSize,I);
%MATLAB_Gamma = trace(logm(A))
[ldUbaru,z1] = Lanc_Quad_LogDet(A,m,nvecs);
[hodlr_ld,z1] = SLQ(@(b) HODLRMatVec(K,b),n,m,nvecs);
[ld,z1] = SLQ(kMtrxFcn,n,m,nvecs);
figure()
hold on
plot(hodlr_ld,"HODLR")
plot(ld,"SLQ")
plot(ldUbaru, "Ubaru")