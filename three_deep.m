% TODO: 
% Add documentation
n = 1000; diag_size = 130;
% k = # of samples; approximately the rank of each small matrix
k = 10; % what should this rank really be? We don't know 
I = [1 n];
% make binary tree
idx_tree = tree(I);
idx_tree = createChildren(idx_tree,1,diag_size);
% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
% Pointervals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
K_mtrx = rbf(ep,DM);

% store U and Z values in nodes respective to index intervals.
U_tree = idx_tree;
Z_tree = idx_tree;
% convert # of elements to # of layers; 

% Layer 1 
% create a random matrix
omega = randn(n,2*k);
% interlace random matrix 
omega(1:500,1:k) = 0; omega(501:end,k+1:end) = 0;
% create sample matrix
Y = K_mtrx * omega;
% index the orthogonalized sample matrix for relevant blocks
U2 = orth(Y(1:500,1:k)); U3 = orth(Y(501:end,k+1:end));
% allocate space for a 'U' matrix
U_mtrx = zeros(n,2*k);
% interlace orthonormal U matrix
U_mtrx(1:500,k+1:end) = U2;
U_mtrx(501:end,1:k) = U3; 
% create the Z matrix
Z_mtrx = K_mtrx' * U_mtrx;
% retrieve the relevant blocks 
Z2 = Z_mtrx(1:500,1:k); Z3 = Z_mtrx(501:end,k+1:end);

% K2
U_tree = U_tree.set(2,U2);
Z_tree = Z_tree.set(2,Z2);

% K3
U_tree = U_tree.set(3,U3);
Z_tree = Z_tree.set(3,Z3);

%disp(abs(norm(K_mtrx(1:500,501:end) - (U_tree.get(2) * Z_tree.get(3)'),'fro') / norm(K_mtrx(1:500,501:end),'fro')));
%disp(abs(norm(K_mtrx(501:end,1:500) - (U_tree.get(3) * Z_tree.get(2)'),'fro') / norm(K_mtrx(501:end,1:500),'fro')));

% Layer 2 (4 blocks)
% create a random matrix
omega = randn(n,2*k);
% interlace random matrix
omega(1:250,1:k) = 0; omega(251:500,k+1:end) = 0;
omega(501:750,1:k) = 0; omega(751:end,k+1:end) = 0;
% create sample matrix
Y = K_mtrx * omega;
% pre-allocate space for upcoming result 
Y4 = Y(1:250,1:k); Y5 = Y(251:500,k+1:end);
Y6 = Y(501:750,1:k); Y7 = Y(751:end,k+1:end);

% regardless of B or not, we get the same result.
%% B: temporary variable; often used to store residual values %%
% subtract away interactions irrelevant to this level.
B = table(U_tree.get(2)).Var1(1:250,:) * (table(Z_tree.get(3)).Var1(251:end,:)' * omega(751:end,1:k));
%B2 = table(U_tree.get(2) * (Z_tree.get(3)' * omega(501:end,1:k))).Var1(1:250,:);
U4 = orth(Y4 - B);

B = table(U_tree.get(2)).Var1(251:end,:) * (table(Z_tree.get(3)).Var1(1:250,:)' * omega(501:750,k+1:end));
%B = table((U_tree.get(2) * Z_tree.get(3)') * omega(501:end,k+1:end)).Var1(251:end,:);
U5 = orth(Y5 - B);

B = table(U_tree.get(3)).Var1(1:250,:) * (table(Z_tree.get(2)).Var1(251:end,:)' * omega(251:500,1:k));
%B = table((U_tree.get(3) * Z_tree.get(2)') * omega(1:500,1:k)).Var1(1:250,:);
U6 = orth(Y6 - B);

B = table(U_tree.get(3)).Var1(251:end,:) * (table(Z_tree.get(2)).Var1(1:250,:)' * omega(1:250,k+1:end));
%B = table((U_tree.get(3) * Z_tree.get(2)') * omega(1:500,k+1:end)).Var1(251:end,:);
U7 = orth(Y7 - B);
% preallocate space for U matrix
U_mtrx = zeros(n,2*k);
% interlace to form U matrix
%% this only works for k = 9 or smaller.
U_mtrx(1:250,k+1:end) = U4; U_mtrx(251:500,1:k) = U5;
U_mtrx(501:750,k+1:end) = U6; U_mtrx(751:end,1:k) = U7;
% create Z matrix
Z_mtrx = K_mtrx' * U_mtrx;
% retrieve only the relevant portions
Z4 = Z_mtrx(1:250,1:k); Z5 = Z_mtrx(251:500,k+1:end);
Z6 = Z_mtrx(501:750,1:k); Z7 = Z_mtrx(751:end,k+1:end);
% (AB^T)^T = BA^T
B = table(Z_tree.get(2)).Var1(1:250,:) * (table(U_tree.get(3)).Var1(251:end,:)' * U_mtrx(751:end,1:k));
%B2 = table(U_tree.get(2) * (Z_tree.get(3)' * omega(501:end,1:k))).Var1(1:250,:);
Z4 = Z4 - B;

B = table(Z_tree.get(2)).Var1(251:end,:) * (table(U_tree.get(3)).Var1(1:250,:)' * U_mtrx(501:750,k+1:end));
Z5 = Z5 - B;
K_exact = U4 * Z5';
size(K_exact)
return

B = table(Z_tree.get(3)).Var1(1:250,:) * (table(U_tree.get(2)).Var1(251:end,:)' * U_mtrx(251:500,1:k));
Z6 = Z6 - B;

B = table(Z_tree.get(3)).Var1(251:end,:) * (table(U_tree.get(2)).Var1(1:250,:)' * U_mtrx(1:250,k+1:end));
Z7 = Z7 - B;
%disp(abs(norm(K_mtrx(1:250,251:500) - (U4 * Z5'), 'fro') / norm(K_mtrx(1:250,251:500),'fro')));
%disp(abs(norm(K_mtrx(251:500,1:250) - (U5 * Z4'), 'fro') / norm(K_mtrx(251:500,1:250),'fro')));
%disp(abs(norm(K_mtrx(501:750,751:end) - (U6 * Z7'), 'fro') / norm(K_mtrx(501:750,751:end),'fro')));
%disp(abs(norm(K_mtrx(751:end,501:750) - (U7 * Z6'), 'fro') / norm(K_mtrx(751:end,501:750),'fro')));

% K4
U_tree = U_tree.set(4,U4);
Z_tree = Z_tree.set(4,Z4);
% K5
U_tree = U_tree.set(5,U5);
Z_tree = Z_tree.set(5,Z5);
% K6
U_tree = U_tree.set(6,U6);
Z_tree = Z_tree.set(6,Z6);
% K7
U_tree = U_tree.set(7,U7);
Z_tree = Z_tree.set(7,Z7);

% Layer 3 (8 blocks)
% create a random matrix
omega = randn(n,2*k);
% interlace random matrix
omega(1:125,1:k) = 0; omega(126:250,k+1:end) = 0;
omega(251:375,1:k) = 0; omega(376:500,k+1:end) = 0;
omega(501:625,1:k) = 0; omega(626:750,k+1:end) = 0;
omega(751:875,1:k) = 0; omega(876:end,k+1:end) = 0; 
% create sample matrix
Y = K_mtrx * omega;
% allocate "impure" Y samples.
U8 = Y(1:125,1:k); U9 = Y(126:250,k+1:end);
U10 = Y(251:375,1:k); U11 = Y(376:500,k+1:end);
U12 = Y(501:625,1:k); U13 = Y(626:750,k+1:end); 
U14 = Y(751:875,1:k); U15 = Y(876:end,k+1:end);

% Purify and orthogonalize Y samples
B = table((U_tree.get(2) * Z_tree.get(3)') * omega(501:end,1:k)).Var1(1:125,:) ...
+ table((U_tree.get(4) * Z_tree.get(5)') * omega(251:500,1:k)).Var1(1:125,:);
U8 = orth(U8 - B); 

B = table((U_tree.get(2) * Z_tree.get(3)') * omega(501:end,k+1:end)).Var1(126:250,:)...
+ table((U_tree.get(4) * Z_tree.get(5)') * omega(251:500,k+1:end)).Var1(126:250,:);
U9 = orth(U9 - B);

B = table((U_tree.get(2) * Z_tree.get(3)') * omega(501:end,1:k)).Var1(251:375,:)...
+ table((U_tree.get(5) * Z_tree.get(4)') * omega(1:250,1:k)).Var1(1:125,:);
U10 = orth(U10 - B);

B = table(( U_tree.get(2) * Z_tree.get(3)') * omega(501:end,k+1:end)).Var1(376:500,:)...
+ table((U_tree.get(5) * Z_tree.get(4)') * omega(1:250,k+1)).Var1(126:250,:);
U11 = orth(U11 - B);

B = table((U_tree.get(3) * Z_tree.get(2)') * omega(1:500,1:k)).Var1(1:125,:)...
+ table((U_tree.get(6) * Z_tree.get(7)') * omega(751:end,1:k)).Var1(1:125,:);
U12 = orth(U12- B); 

B = table((U_tree.get(3) * Z_tree.get(2)') * omega(1:500,k+1:end)).Var1(126:250,:)...
+ table((U_tree.get(6) * Z_tree.get(7)') * omega(751:end,k+1:end)).Var1(126:250,:);
U13 = orth(U13 - B);

B = table((U_tree.get(3) * Z_tree.get(2)') * omega(1:500,1:k)).Var1(251:375,:)...
+ table(( U_tree.get(7) * Z_tree.get(6)') * omega(501:750,1:k)).Var1(1:125,:);
U14 = orth(U14 - B);

B = table((U_tree.get(3) * Z_tree.get(2)') * omega(1:500,k+1:end)).Var1(376:500,:)...
+ table((U_tree.get(7) *Z_tree.get(6)') * omega(501:750,k+1:end)).Var1(126:250,:);
U15 = orth(U15 - B);

% preallocate for U matrix
U_mtrx = zeros(n,2*k);
% interlace to form U matrix
U_mtrx(1:125,k+1:end) = U8; U_mtrx(126:250,1:k) = U9;
U_mtrx(251:375,k+1:end) = U10; U_mtrx(376:500,1:k) = U11;
U_mtrx(501:625,k+1:end) = U12; U_mtrx(626:750,1:k) = U13;
U_mtrx(751:875,k+1:end) = U14; U_mtrx(876:end,1:k) = U15;
% create Z matrix
Z_mtrx = K_mtrx' * U_mtrx;

% retrieve the relevant blocks
Z8 = Z_mtrx(1:125,1:k); Z9 = Z_mtrx(126:250,k+1:end);
Z10 = Z_mtrx(251:375,1:k); Z11 = Z_mtrx(376:500,k+1:end);
Z12 = Z_mtrx(501:625,1:k); Z13 = Z_mtrx(626:750,k+1:end);
Z14 = Z_mtrx(751:875,1:k); Z15 = Z_mtrx(876:end,k+1:end);

B = table((U_tree.get(3) * Z_tree.get(2)')' * U_mtrx(501:end,1:k)).Var1(1:125,:) ...
+ table((U_tree.get(5) * Z_tree.get(4)')' * U_mtrx(251:500,1:k)).Var1(1:125,:);
Z8 = Z8 - B; 

B = table((U_tree.get(3) * Z_tree.get(2)')' * U_mtrx(501:end,k+1:end)).Var1(126:250,:)...
+ table((U_tree.get(5) * Z_tree.get(4)')' * U_mtrx(251:500,k+1:end)).Var1(126:250,:);
Z9 = Z9 - B;

B = table((U_tree.get(3) * Z_tree.get(2)')' * U_mtrx(501:end,1:k)).Var1(251:375,:)...
+ table((U_tree.get(4) * Z_tree.get(5)')' * U_mtrx(1:250,1:k)).Var1(1:125,:);
Z10 = Z10 - B;

B = table((U_tree.get(3) * Z_tree.get(2)')' * U_mtrx(501:end,k+1:end)).Var1(376:500,:)...
+ table((U_tree.get(4) * Z_tree.get(5)')' * U_mtrx(1:250,k+1)).Var1(126:250,:);
Z11 = Z11 - B;

B = table((U_tree.get(2) * Z_tree.get(3)')' * U_mtrx(1:500,1:k)).Var1(1:125,:)...
+ table((U_tree.get(7) * Z_tree.get(6)')' * U_mtrx(751:end,1:k)).Var1(1:125,:);
Z12 = Z12- B;

B = table((U_tree.get(2) * Z_tree.get(3)')' * U_mtrx(1:500,k+1:end)).Var1(126:250,:)...
+ table((U_tree.get(7) * Z_tree.get(6)')' * U_mtrx(751:end,k+1:end)).Var1(126:250,:);
Z13 = Z13 - B;

B = table((U_tree.get(2) * Z_tree.get(3)')' * U_mtrx(1:500,1:k)).Var1(251:375,:)...
+ table((U_tree.get(6) * Z_tree.get(7)')' * U_mtrx(501:750,1:k)).Var1(1:125,:);
Z14 = Z14 - B;

B = table(( U_tree.get(2) * Z_tree.get(3)')' * U_mtrx(1:500,k+1:end)).Var1(376:500,:)...
+ table((U_tree.get(6) *Z_tree.get(7)')' * U_mtrx(501:750,k+1:end)).Var1(126:250,:);
Z15 = Z15 - B;

%disp(abs(norm(K_mtrx(1:125,126:250) - (U8 * Z9'), 'fro')...
%/ norm(K_mtrx(1:125,126:250),'fro')));
%disp(abs(norm(K_mtrx(126:250,1:125) - (U9 * Z8'), 'fro')...
%/ norm(K_mtrx(126:250,1:125),'fro')));
%disp(abs(norm(K_mtrx(251:375,376:500) - (U10 * Z11'), 'fro')...
%/ norm(K_mtrx(251:375,375:500),'fro')));
%disp(abs(norm(K_mtrx(376:500,251:375) - (U11 * Z10'), 'fro')...
%/ norm(K_mtrx(376:500,251:375),'fro')));
%disp(abs(norm(K_mtrx(501:625,626:750) - (U12 * Z13'), 'fro')...
%/ norm(K_mtrx(501:625,626:750),'fro')));
%disp(abs(norm(K_mtrx(626:750,501:625) - (U13 * Z12'), 'fro')...
%/ norm(K_mtrx(626:750,501:625),'fro')));
%disp(abs(norm(K_mtrx(751:875,876:end) - (U14 * Z15'), 'fro')...
%/ norm(K_mtrx(751:875,875:end),'fro')));
%disp(abs(norm(K_mtrx(876:end,751:875) - (U15 * Z14'), 'fro')...
%/ norm(K_mtrx(876:end,751:875),'fro')));
% K8
U_tree = U_tree.set(8,U8);
Z_tree = Z_tree.set(8,Z8);

% K9
U_tree = U_tree.set(9,U9);
Z_tree = Z_tree.set(9,Z9);

% K10
U_tree = U_tree.set(10,U10);
Z_tree = Z_tree.set(10,Z10);

% K11
U_tree = U_tree.set(11,U11);
Z_tree = Z_tree.set(11,Z11);

% K12
U_tree = U_tree.set(12,U12);
Z_tree = Z_tree.set(12,Z12);

% K13
U_tree = U_tree.set(13,U13);
Z_tree = Z_tree.set(13,Z13);

% K14
U_tree = U_tree.set(14,U14);
Z_tree = Z_tree.set(14,Z14);

% K15
U_tree = U_tree.set(15,U15);
Z_tree = Z_tree.set(15,Z15);

% reconstruct the K matrix approximation from U and Z values
K_approx = K_mtrx;
K_approx(1:500,501:end) = U_tree.get(2) * Z_tree.get(3)';
K_approx(501:end,1:500) = U_tree.get(3) * Z_tree.get(2)';

K_approx(1:250,251:500) = U_tree.get(4) * Z_tree.get(5)';
K_approx(251:500,1:250) = U_tree.get(5) * Z_tree.get(4)';
K_approx(501:750,751:end) = U_tree.get(6) * Z_tree.get(7)';
K_approx(751:end,501:750) = U_tree.get(7) * Z_tree.get(6)';

K_approx(1:125,126:250) = U_tree.get(8) * Z_tree.get(9)';
K_approx(126:250,1:125) = U_tree.get(9) * Z_tree.get(8)';
K_approx(251:375,376:500) = U_tree.get(10) * Z_tree.get(11)';
K_approx(376:500,251:375) = U_tree.get(11) * Z_tree.get(10)';
K_approx(501:625,626:750) = U_tree.get(12) * Z_tree.get(13)';
K_approx(626:750,501:625) = U_tree.get(13) * Z_tree.get(12)';
K_approx(751:875,876:end) = U_tree.get(14) * Z_tree.get(15)';
K_approx(876:end,751:875) = U_tree.get(15) * Z_tree.get(14)';

%imagesc(K_mtrx - K_approx)
disp(abs(norm(K_mtrx - K_approx, 'fro') / norm(K_mtrx, 'fro')))
%end
