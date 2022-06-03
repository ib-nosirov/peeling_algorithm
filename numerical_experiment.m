% TODO: 
% Make it so we only multiply by K.
% Add documentation
n = 1000;
diag_size = 130;
I = [1 n];
% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 5; dim = 1;
% Pointervals at which to sample
x = CreatePoints(n^dim,dim,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
K_mtrx = rbf(ep,DM);

% k = # of samples; approximately the rank of each small matrix
k = 10; % what should this rank really be? We don't know 

% store U and Z values in nodes respective to index intervals.
% make binary tree
idx_tree = tree(I);
idx_tree = create_children(idx_tree,1,diag_size);
t_order = idx_tree.breadthfirstiterator
U_tree = idx_tree;
Z_tree = idx_tree;

% convert # of elements to # of layers; 
%ex: 15 nodes -> 3 layers
tree_depth = log2(nnodes(idx_tree)+1) - 1;
Omega = randn(n, 2*k);
% Layer 1 (2 blocks)

for l=2:tree_depth
    % create a random matrix
    omega = randn(n,2*k);
    % interlace random matrix 
    mask = Interlace(idx_tree,l,k);
    % compute sample matrix using interlaced random matrix.
    Y = K_mtrx * (Omega .* mask);
    % partition size -> # of intervals
    part_size = floor(n/(2^l)); intervals = 1:part_size:n+1;
    % orthogonalize and separate into skinny matrices
    for idx = 1:length(intervals)-1
        % upper diagonal submatrices always have even threshold values
        quotient = floor(intervals(idx) / part_size);
        % top diagonal K submatrices
        if mod(quotient,2) == 0
            % orthonormalize each green block in Fig 12a and place it in U_tree
            Y_pure = PurifySample(Y,intervals,idx); 
	    U_tree = U_tree.set(t_order(it),orth(Y_pure));
            % Perform the operation in Fig 13 and store the result in Z_tree
            % multiply by the transpose of K_23
	    Z_tree = Z_tree.set(t_order(it), K' * U_tree.get(t_order(it)));
            % zero out submatrices in K
        % bottom diagonal K submatrices
        else
            % orthonormalize each green block in Fig 12a and place it in U_tree
            U_tree = U_tree.set(t_order(it),...
                     orth(Y(intervals(idx):intervals(idx)+part_size-1, k+1:end)));
            % Perform the operation in Fig 13 and store the result in Z_tree
            Z_tree = Z_tree.set(t_order(it), K' * U_tree.get(t_order(it)));
        end
        it = it + 1;
    end

%    function Y_pure = PurifySample(Y,intervals,idx)
%        
%    Y(intervals(idx):intervals(idx+1)-1, 1:k);
%    end
end
%% check the accuracy.
%K = rbf(ep, DM);
%for i = 3:2:length(idx_tree.breadthfirstiterator)
%  idx = t_order(i-1);
%  idx_next = t_order(i);
%  % upper triangle
%  x_idx = idx_tree.get(idx); y_idx = idx_tree.get(idx_next);
%  x_idx
%  y_idx
%  true_K = K(x_idx(1):x_idx(2), y_idx(1):y_idx(2));
%  % MSE
%  disp(norm(true_K - U_tree.get(idx) * (Z_tree.get(idx))', 'fro')^2 / numel(true_K) / norm(U_tree.get(idx) * (Z_tree.get(idx))',"fro")^2)
%  % lower triangle
%  x_idx = idx_tree.get(idx_next); y_idx = idx_tree.get(idx);
%  x_idx
%  y_idx
%  true_K = K(x_idx(1):x_idx(2), y_idx(1):y_idx(2));
%  % MSE
%  disp(norm(true_K - U_tree.get(idx_next) * (Z_tree.get(idx_next))', 'fro')^2 / numel(true_K) / norm(U_tree.get(idx_next) * (Z_tree.get(idx_next))',"fro")^2)
%end

function mask = MakeMask(n,l,k)
    % build the first "row" of zeros and random values (Blue/white block in
    % Fig (12a)
    row = repmat([zeros(1,k) ones(1,k)], floor(n/(2^l)), 1);
    % Flip the existing block and stack the two to form a "block"
    row_flip = fliplr(row);
    block = [row; row_flip];
    % keep stacking blocks 2^(l-1) times to form the entire mask
    mask = repmat(block, 2^(l-1), 1);
end

function interlaced_matrix = Interlace(tree,l,k)
end
