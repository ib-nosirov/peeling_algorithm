% TODO: 
% Make it so we only multiply by K.
% Add documentation
n = 1000;
d = 130;
% k = # of samples; approximately the rank of each small matrix
k = 10; % what should this rank really be? We don't know 
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
t_order = idx_tree.breadthfirstiterator
it = 2;
U_tree = idx_tree;
Z_tree = idx_tree;
% convert # of elements to # of layers; 

% Layer 1 
omega = randn(n,2*k);
omega(1:500,1:k) = 0; omega(501:1000,k+1:end) = 0;
Y = K_mtrx * omega;
U = orth(Y);
U_mtrx = zeros(n,2*k);
U_top = U(1:500,1:k);
U_bottom = U(501:end,k+1:end);
U_mtrx(1:500,1:k) = U_top; U_mtrx(501:end,k+1:end) = U_bottom; 
Z_mtrx = K_mtrx' * U_mtrx;
Z_top = Z_mtrx(1:500,1:k);
Z_bottom = Z_mtrx(501:end,k+1:end);

% K23
U_tree = U_tree.set(t_order(2),U_top);
Z_tree = Z_tree.set(t_order(2),Z_top);

% K32
U_tree = U_tree.set(t_order(3),U_bottom);
Z_tree = Z_tree.set(t_order(3),Z_bottom);

% Layer 2 
omega = randn(n,2*k);
omega(1:250,1:k) = 0; omega(251:500,k+1:end) = 0; omega(501:750,1:k) = 0;
omega(751:1000,k+1:end) = 0;
Y = K_mtrx * omega;
U45 = Y(1:250,1:k); U54 = Y(251:500,k+1:end); U67 = Y(501:750,1:k);
U76 = Y(751:1000,k+1:end);

tmp = table(( U_tree.get(2) * Z_tree.get(2)' ) *...
omega(501:1000,1:k)).Var1(1:250,:);
U45 = orth(U45 - tmp); 
tmp = table(( U_tree.get(2) * Z_tree.get(2)' ) *...
omega(501:1000,k+1:end)).Var1(251:end,:);
U54 = orth(U54 - tmp);
tmp = table(( U_tree.get(3) * Z_tree.get(3)' ) *...
omega(1:500,1:k)).Var1(1:250,:);
U67 = orth(U67 - tmp);
tmp = table(( U_tree.get(3) * Z_tree.get(3)' ) *...
omega(1:500,k+1:end)).Var1(251:end,:);
U76 = orth(U76 - tmp);

U_mtrx = zeros(n,2*k);
U_mtrx(1:250,1:k) = U45; U_mtrx(251:500,k+1:end) = U54;
U_mtrx(501:750,1:k) = U67; U_mtrx(751:1000,k+1:end) = U76;

Z_mtrx = K_mtrx' * U_mtrx;
Z45 = Z_mtrx(1:250,1:k); Z54 = Z_mtrx(251:500,k+1:end);
Z67 = Z_mtrx(501:750,1:k); Z76 = Z_mtrx(751:1000,k+1:end);

% K45
U_tree = U_tree.set(t_order(4),U45);
Z_tree = Z_tree.set(t_order(4),Z45);

% K54
U_tree = U_tree.set(t_order(5),U54);
Z_tree = Z_tree.set(t_order(5),Z54);

% K67
U_tree = U_tree.set(t_order(6),U67);
Z_tree = Z_tree.set(t_order(6),Z67);

% K76
U_tree = U_tree.set(t_order(7),U76);
Z_tree = Z_tree.set(t_order(7),Z76);

% Layer 3
omega = randn(n,2*k);
omega(1:125,1:k) = 0; omega(126:250,k+1:end) = 0; omega(251:375,1:k) = 0;
omega(376:500,k+1:end) = 0; omega(501:625,1:k) = 0; omega(626:750,k+1:end) = 0;
omega(751:875,1:k) = 0; omega(876:1000,k+1:end) = 0; 

Y = K_mtrx * omega;
U89 = Y(1:125,1:k); U98 = Y(126:250,k+1:end); U1011 = Y(251:375,1:k);
U1110 = Y(376:500,k+1:end); U1213 = Y(501:625,1:k); U1312 = Y(626:750,k+1:end); 
U1415 = Y(751:875,1:k); U1514 = Y(876:1000,k+1:end);

tmp = table(( U_tree.get(2) * Z_tree.get(2)' ) *...
omega(501:1000,1:k)).Var1(1:125,:) + table((U_tree.get(4) * Z_tree.get(4)')...
* omega(251:500,1:k)).Var1(1:125,:);
U89 = orth(U89 - tmp); 
tmp = table(( U_tree.get(2) * Z_tree.get(2)' ) *...
omega(501:1000,k+1:end)).Var1(126:250,:) + table((U_tree.get(4) * Z_tree.get(4)') *...
omega(251:500,k+1:end)).Var1(126:250,:);
U98 = orth(U98 - tmp);
tmp = table(( U_tree.get(2) * Z_tree.get(2)' ) *...
omega(501:1000,1:k)).Var1(251:375,:) + table((U_tree.get(5) * Z_tree.get(5)') *...
omega(1:250,1:k)).Var1(1:125,:);
U1011 = orth(U1011 - tmp);
tmp = table(( U_tree.get(2) * Z_tree.get(2)' ) *...
omega(501:1000,k+1:end)).Var1(376:500,:) + table((U_tree.get(5) * Z_tree.get(5)')...
* omega(1:250,k+1)).Var1(126:250,:);
U1110 = orth(U1110 - tmp);

% complete the tmp
tmp = table(( U_tree.get(3) * Z_tree.get(3)' ) *...
omega(1:500,1:k)).Var1(1:125,:) + table(( U_tree.get(6) * Z_tree.get(6)' )...
* omega(751:1000,1:k)).Var1(1:125,:);
U1213 = orth(U1213- tmp); 
tmp = table((U_tree.get(3) * Z_tree.get(3)') *...
omega(1:500,k+1:end)).Var1(126:250,:) + table(( U_tree.get(6) * Z_tree.get(6)' ) *...
omega(751:1000,k+1:end)).Var1(126:250,:);
U1312 = orth(U1312 - tmp);
tmp = table((U_tree.get(3) * Z_tree.get(3)') *...
omega(1:500,1:k)).Var1(251:375,:) + table(( U_tree.get(7) * Z_tree.get(7)')...
* omega(501:750,1:k)).Var1(1:125,:);
U1415 = orth(U1415 - tmp);
tmp = table(( U_tree.get(3) * Z_tree.get(3)' ) *...
omega(1:500,k+1:end)).Var1(376:500,:) + table((U_tree.get(7) *Z_tree.get(7)')...
* omega(501:750,k+1:end)).Var1(126:250,:);
U1514 = orth(U1514 - tmp);

U_mtrx = zeros(n,2*k);
U_mtrx(1:125,k+1:end) = U89;
U_mtrx(126:250,1:k) = U98;
U_mtrx(251:375,k+1:end) = U1011;
U_mtrx(376:500,1:k) = U1110;
U_mtrx(501:625,k+1:end) = U1213;
U_mtrx(626:750,1:k) = U1312;
U_mtrx(751:875,k+1:end) = U1415;
U_mtrx(876:1000,1:k) = U1514;

Z_mtrx = K_mtrx' * U_mtrx;

Z89 = Z_mtrx(1:125,1:k);
Z98 = Z_mtrx(126:250,k+1:end);
Z1011 = Z_mtrx(251:375,1:k);
Z1110 = Z_mtrx(376:500,k+1:end);
Z1213 = Z_mtrx(501:625,1:k);
Z1312 = Z_mtrx(626:750,k+1:end);
Z1415 = Z_mtrx(751:875,1:k);
Z1514 = Z_mtrx(876:1000,k+1:end);

% K89
U_tree = U_tree.set(t_order(8),U89);
Z_tree = Z_tree.set(t_order(8),Z89);

% K98
U_tree = U_tree.set(t_order(9),U98);
Z_tree = Z_tree.set(t_order(9),Z98);

% K1011
U_tree = U_tree.set(t_order(10),U1011);
Z_tree = Z_tree.set(t_order(10),Z1011);

% K1110
U_tree = U_tree.set(t_order(11),U1110);
Z_tree = Z_tree.set(t_order(11),U1110);

% K1213
U_tree = U_tree.set(t_order(12),U1213);
Z_tree = Z_tree.set(t_order(12),Z1213);

% K1312
U_tree = U_tree.set(t_order(13),U1312);
Z_tree = Z_tree.set(t_order(13),Z1312);

% K1415
U_tree = U_tree.set(t_order(14),U1415);
Z_tree = Z_tree.set(t_order(14),Z1415);

% K1514
U_tree = U_tree.set(t_order(15),U1514);
Z_tree = Z_tree.set(t_order(15),Z1514);

K_approx = K_mtrx;
K_approx(1:500,501:1000) = U_tree.get(2) * Z_tree.get(2)';
K_approx(501:1000,1:500) = U_tree.get(3) * Z_tree.get(3)';
K_approx(1:250,251:500) = U_tree.get(4) * Z_tree.get(4)';
K_approx(251:500,1:250) = U_tree.get(5) * Z_tree.get(5)';
K_approx(501:750,751:1000) = U_tree.get(10) * Z_tree.get(10)';
K_approx(751:1000,501:750) = U_tree.get(11) * Z_tree.get(11)';
K_approx(1:125,126:250) = U_tree.get(6) * Z_tree.get(6)';
K_approx(126:250,1:125) = U_tree.get(7) * Z_tree.get(7)';
K_approx(251:375,376:500) = U_tree.get(8) * Z_tree.get(8)';
K_approx(376:500,251:375) = U_tree.get(9) * Z_tree.get(9)';
K_approx(501:625,626:750) = U_tree.get(12) * Z_tree.get(12)';
K_approx(626:750,501:625) = U_tree.get(13) * Z_tree.get(13)';
K_approx(751:875,876:1000) = U_tree.get(14) * Z_tree.get(14)';
K_approx(876:1000,751:875) = U_tree.get(15) * Z_tree.get(15)';


disp(norm(K_mtrx - K_approx, 'fro')^2 / numel(K_mtrx))
%end
