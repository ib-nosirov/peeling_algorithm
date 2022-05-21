n = 100;
diag_size = 13;
I = [1 n];
% make binary tree.
idx_tree = tree(I);
idx_tree = create_children(idx_tree, 1, diag_size);
%disp(idx_tree.tostring);
% make a kernel.
rbf = @(e,r) exp(-(e*r).^2); ep = 5; d = 1;
% Points at which to sample.
x = CreatePoints(n^d, d, 'u');
%x = x(randperm(length(x)));
% compute the absolute difference.
DM = DistanceMatrix(x, x);
% sample matrix.
K = rbf(ep, DM);

% k = # of samples; approximately the rank of each small matrix
k = 10; % what should this rank really be? We don't know 

% store U and Z values in nodes respective to index intervals.
t_order = idx_tree.breadthfirstiterator;
it = 2;
U_tree = idx_tree;
Z_tree = idx_tree;

for l=1:3
  % Omega = random matrix.
  Omega = randn(n, 2*k);
  % l = generation in the tree
  thd = floor(n/(2^l));
  % make intervals (ints)
  ints = 1:thd:n+1;
  % build mask
  mask = 1:n;
  for idx = 1:length(ints)-2
    mask(ints(idx):ints(idx+1)-1) = mod(mask(ints(idx):ints(idx+1)-1) / thd, 2) == 0;
  end
% there is something wrong here.
  mask = mask ~= 0; mask = [repmat(mask, k, 1); repmat(flip(mask), k, 1)];
  % interlace
  Omega = Omega .* mask';
  % compute the sample matrix Y.
  Y = K * Omega;
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
