% TODO: 
% Make it so we only multiply by K.
% Add documentation
n = 1000;
d = 130;
% k = # of samples; approximately the rank of each small matrix
k = 9; % what should this rank really be? We don't know 
I = [1 n];
% make binary tree
idx_tree = tree(I);
idx_tree = create_children(idx_tree,1,d);
% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 5; d = 1;
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
U23 = orth(Y(1:500,1:k)); U32 = orth(Y(501:end,k+1:end));
% allocate space for a 'U' matrix
U_mtrx = zeros(n,2*k);
% interlace orthonormal U matrix
U_mtrx(1:500,k+1:end) = U23; U_mtrx(501:end,1:k) = U32; 
% create the Z matrix
Z_mtrx = K_mtrx' * U_mtrx;
% retrieve the relevant blocks 
Z23 = Z_mtrx(1:500,1:k); Z32 = Z_mtrx(501:end,k+1:end);

disp(norm(K_mtrx(1:500,501:end) - (U23 * Z23'), 'fro')^2 / numel(K_mtrx(1:500,501:end)))
return
% K23
U_tree = U_tree.set(2,U23);
Z_tree = Z_tree.set(2,Z23);

% K32
U_tree = U_tree.set(3,U32);
Z_tree = Z_tree.set(3,Z32);

% Layer 2 (4 blocks)
% create a random matrix
omega = randn(n,2*k);
% interlace random matrix
omega(1:250,1:k) = 0; omega(251:500,k+1:end) = 0;
omega(501:750,1:k) = 0; omega(751:1000,k+1:end) = 0;
% create sample matrix
Y = K_mtrx * omega;
% orthogonalize the result
U45 = Y(1:250,1:k); U54 = Y(251:500,k+1:end);
U67 = Y(501:750,1:k); U76 = Y(751:1000,k+1:end);

%% B: temporary variable; often used to store 'B' values %%
% subtract away the irrelevant interactions at previous levels
B = table(( U_tree.get(2) * Z_tree.get(2)' ) * omega(501:1000,1:k)).Var1(1:250,:);
U45 = orth(U45 - B); 

B = table(( U_tree.get(2) * Z_tree.get(2)' ) * omega(501:1000,k+1:end)).Var1(251:end,:);
U54 = orth(U54 - B);

B = table(( U_tree.get(3) * Z_tree.get(3)' ) * omega(1:500,1:k)).Var1(1:250,:);
U67 = orth(U67 - B);

B = table(( U_tree.get(3) * Z_tree.get(3)' ) * omega(1:500,k+1:end)).Var1(251:end,:);
U76 = orth(U76 - B);

% preallocate space for U matrix
U_mtrx = zeros(n,2*k);
% interlace to form U matrix
U_mtrx(1:250,1:k) = U45; U_mtrx(251:500,k+1:end) = U54;
U_mtrx(501:750,1:k) = U67; U_mtrx(751:1000,k+1:end) = U76;
% create Z matrix
Z_mtrx = K_mtrx' * U_mtrx;
% retrieve only the relevant portions
Z45 = Z_mtrx(1:250,1:k); Z54 = Z_mtrx(251:500,k+1:end);
Z67 = Z_mtrx(501:750,1:k); Z76 = Z_mtrx(751:1000,k+1:end);

% K45
U_tree = U_tree.set(4,U45);
Z_tree = Z_tree.set(4,Z45);
% K54
U_tree = U_tree.set(5,U54);
Z_tree = Z_tree.set(5,Z54);
% K67
U_tree = U_tree.set(6,U67);
Z_tree = Z_tree.set(6,Z67);
% K76
U_tree = U_tree.set(7,U76);
Z_tree = Z_tree.set(7,Z76);

% Layer 3 (8 blocks)
% create a random matrix
omega = randn(n,2*k);
% interlace random matrix
omega(1:125,1:k) = 0; omega(126:250,k+1:end) = 0;
omega(251:375,1:k) = 0; omega(376:500,k+1:end) = 0;
omega(501:625,1:k) = 0; omega(626:750,k+1:end) = 0;
omega(751:875,1:k) = 0; omega(876:1000,k+1:end) = 0; 
% create sample matrix
Y = K_mtrx * omega;
% allocate "impure" Y samples.
U89 = Y(1:125,1:k); U98 = Y(126:250,k+1:end);
U1011 = Y(251:375,1:k); U1110 = Y(376:500,k+1:end);
U1213 = Y(501:625,1:k); U1312 = Y(626:750,k+1:end); 
U1415 = Y(751:875,1:k); U1514 = Y(876:1000,k+1:end);

% Purity and orthogonalize Y samples
B = table(( U_tree.get(2) * Z_tree.get(2)' ) * omega(501:1000,1:k)).Var1(1:125,:) ...
+ table((U_tree.get(4) * Z_tree.get(4)') * omega(251:500,1:k)).Var1(1:125,:);
U89 = orth(U89 - B); 

B = table(( U_tree.get(2) * Z_tree.get(2)' ) * omega(501:1000,k+1:end)).Var1(126:250,:)...
+ table((U_tree.get(4) * Z_tree.get(4)') * omega(251:500,k+1:end)).Var1(126:250,:);
U98 = orth(U98 - B);

B = table(( U_tree.get(2) * Z_tree.get(2)' ) * omega(501:1000,1:k)).Var1(251:375,:)...
+ table((U_tree.get(5) * Z_tree.get(5)') * omega(1:250,1:k)).Var1(1:125,:);
U1011 = orth(U1011 - B);

B = table(( U_tree.get(2) * Z_tree.get(2)' ) * omega(501:1000,k+1:end)).Var1(376:500,:)...
+ table((U_tree.get(5) * Z_tree.get(5)') * omega(1:250,k+1)).Var1(126:250,:);
U1110 = orth(U1110 - B);

B = table(( U_tree.get(3) * Z_tree.get(3)' ) * omega(1:500,1:k)).Var1(1:125,:)...
+ table(( U_tree.get(6) * Z_tree.get(6)' ) * omega(751:1000,1:k)).Var1(1:125,:);
U1213 = orth(U1213- B); 

B = table((U_tree.get(3) * Z_tree.get(3)') *omega(1:500,k+1:end)).Var1(126:250,:)...
+ table(( U_tree.get(6) * Z_tree.get(6)' ) * omega(751:1000,k+1:end)).Var1(126:250,:);
U1312 = orth(U1312 - B);

B = table((U_tree.get(3) * Z_tree.get(3)') * omega(1:500,1:k)).Var1(251:375,:)...
+ table(( U_tree.get(7) * Z_tree.get(7)') * omega(501:750,1:k)).Var1(1:125,:);
U1415 = orth(U1415 - B);

B = table(( U_tree.get(3) * Z_tree.get(3)' ) * omega(1:500,k+1:end)).Var1(376:500,:)...
+ table((U_tree.get(7) *Z_tree.get(7)') * omega(501:750,k+1:end)).Var1(126:250,:);
U1514 = orth(U1514 - B);

% preallocate for U matrix
U_mtrx = zeros(n,2*k);
% interlace to form U matrix
U_mtrx(1:125,k+1:end) = U89; U_mtrx(126:250,1:k) = U98;
U_mtrx(251:375,k+1:end) = U1011; U_mtrx(376:500,1:k) = U1110;
U_mtrx(501:625,k+1:end) = U1213; U_mtrx(626:750,1:k) = U1312;
U_mtrx(751:875,k+1:end) = U1415; U_mtrx(876:1000,1:k) = U1514;

% create Z matrix
Z_mtrx = K_mtrx' * U_mtrx;

% retrieve the relevant blocks
Z89 = Z_mtrx(1:125,1:k); Z98 = Z_mtrx(126:250,k+1:end);
Z1011 = Z_mtrx(251:375,1:k); Z1110 = Z_mtrx(376:500,k+1:end);
Z1213 = Z_mtrx(501:625,1:k); Z1312 = Z_mtrx(626:750,k+1:end);
Z1415 = Z_mtrx(751:875,1:k); Z1514 = Z_mtrx(876:1000,k+1:end);

% K89
U_tree = U_tree.set(8,U89);
Z_tree = Z_tree.set(8,Z89);

% K98
U_tree = U_tree.set(9,U98);
Z_tree = Z_tree.set(9,Z98);

% K1011
U_tree = U_tree.set(10,U1011);
Z_tree = Z_tree.set(10,Z1011);

% K1110
U_tree = U_tree.set(11,U1110);
Z_tree = Z_tree.set(11,U1110);

% K1213
U_tree = U_tree.set(12,U1213);
Z_tree = Z_tree.set(12,Z1213);

% K1312
U_tree = U_tree.set(13,U1312);
Z_tree = Z_tree.set(13,Z1312);

% K1415
U_tree = U_tree.set(14,U1415);
Z_tree = Z_tree.set(14,Z1415);

% K1514
U_tree = U_tree.set(15,U1514);
Z_tree = Z_tree.set(15,Z1514);

% reconstruct the K matrix approximation from U and Z values
K_approx = K_mtrx;
K_approx(1:500,501:1000) = U_tree.get(2) * Z_tree.get(2)';
K_approx(501:1000,1:500) = U_tree.get(3) * Z_tree.get(3)';
K_approx(1:250,251:500) = U_tree.get(4) * Z_tree.get(4)';
K_approx(251:500,1:250) = U_tree.get(5) * Z_tree.get(5)';
K_approx(501:750,751:1000) = U_tree.get(6) * Z_tree.get(6)';
K_approx(751:1000,501:750) = U_tree.get(7) * Z_tree.get(7)';
K_approx(1:125,126:250) = U_tree.get(8) * Z_tree.get(8)';
K_approx(126:250,1:125) = U_tree.get(9) * Z_tree.get(9)';
K_approx(251:375,376:500) = U_tree.get(10) * Z_tree.get(10)';
K_approx(376:500,251:375) = U_tree.get(11) * Z_tree.get(11)';
K_approx(501:625,626:750) = U_tree.get(12) * Z_tree.get(12)';
K_approx(626:750,501:625) = U_tree.get(13) * Z_tree.get(13)';
K_approx(751:875,876:1000) = U_tree.get(14) * Z_tree.get(14)';
K_approx(876:1000,751:875) = U_tree.get(15) * Z_tree.get(15)';

disp(norm(K_mtrx - K_approx, 'fro')^2 / numel(K_mtrx))
%end
