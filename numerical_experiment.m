n = 1000;
d = 130;
I = [1 n];
% make binary tree
idx_tree = tree(I);
idx_tree = create_children(idx_tree,1,d);
% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 5; d = 1;
% Points at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
K = rbf(ep,DM);

% k = # of samples; approximately the rank of each small matrix
k = 10; % what should this rank really be? We don't know 

% store U and Z values in nodes respective to index intervals.
t_order = idx_tree.breadthfirstiterator;
it = 2;
U_tree = idx_tree;
Z_tree = idx_tree;
tree_depth = log2(nnodes(idx_tree)+1) - 1;
Omega = randn(n, 2*k);

for l=1:tree_depth
  % l = generation in the tree
  thd = floor(n/(2^l));
  % make intervals (ints)
  ints = 1:thd:n+1;
  row = repmat([zeros(1,k) ones(1,k)], floor(n/(2^l)), 1);
  row_flip = fliplr(row);
  block = [row; row_flip];
  mask = repmat(block, 2^(l-1), 1);
  % compute sample matrix using interlaced random matrix.
  Y = K * (Omega .* mask);
  % orthogonalize and separate into skinny matrices
  % populate tree with U and Z vectors.
  for idx = 1:length(ints)-1
    % upper diagonal submatrices always have even threshold values
    quotient = floor(ints(idx) / thd);
    % top diagonal K submatrices
    if mod(quotient, 2) == 0
      U_tree = U_tree.set(t_order(it), orth(Y(ints(idx):ints(idx+1)-1, 1:k)));
      Z_tree = Z_tree.set(t_order(it), K(ints(idx):ints(idx+1)-1, ints(idx+1):ints(idx+2)-1)' * U_tree.get(t_order(it)));
      % zero out submatrices in K
      K(ints(idx):ints(idx+1)-1, ints(idx+1):ints(idx+2)-1) = 0;
    % bottom diagonal K submatrices
    else
      U_tree = U_tree.set(t_order(it), orth(Y(ints(idx):ints(idx)+thd-1, k+1:end)));
      Z_tree = Z_tree.set(t_order(it), K(ints(idx):ints(idx+1)-1, ints(idx-1):ints(idx)-1)' * U_tree.get(t_order(it)));
      % zero out submatrices in K
      K(ints(idx):ints(idx+1)-1, ints(idx-1):ints(idx)-1) = 0;
    end
    it = it + 1;
  end
end
% fix the idx mismatch issue.
%should you store the same U and Z at the same location in their
%respective trees?
K = rbf(ep, DM);
for i = 3:2:length(idx_tree.breadthfirstiterator)
  idx = t_order(i-1);
  idx_next = t_order(i);
  % upper triangle
  x_idx = idx_tree.get(idx); y_idx = idx_tree.get(idx_next);
  x_idx
  y_idx
  true = K(x_idx(1):x_idx(2), y_idx(1):y_idx(2));
  disp(norm(true - U_tree.get(idx) * (Z_tree.get(idx))', 'fro')^2 / numel(true))
  % lower triangle
  x_idx = idx_tree.get(idx_next); y_idx = idx_tree.get(idx);
  x_idx
  y_idx
  true = K(x_idx(1):x_idx(2), y_idx(1):y_idx(2));
  disp(norm(true - U_tree.get(idx_next) * (Z_tree.get(idx_next))', 'fro')^2 / numel(true))
end
