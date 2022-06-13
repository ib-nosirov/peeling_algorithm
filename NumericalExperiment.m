function [U_tree Z_tree idx_tree] = NumericalExperiment(K_mtrx,n,k,diag_size,I)
% make binary tree
idx_tree = tree(I);
idx_tree = CreateChildren(idx_tree,1,diag_size);
% store U and Z values in nodes respective to index intervals.
U_tree = idx_tree;
Z_tree = idx_tree;

% Layer 1
ComputeFirstLayer();

% every layer after layer 1 requires computing residuals from previous layers.
tree_depth = log2(nnodes(idx_tree)+1) - 1;
% subsequent layers
for l = 2:tree_depth
	ComputeNextLayer(l);
end

	function output = ComputeFirstLayer()
	% create a random matrix and interlace it with zeroes.
	omega = randn(n,2*k); omega(1:n/2,1:k) = 0; omega((n/2)+1:end,k+1:end) = 0;
	% create sample matrix
	Y = K_mtrx * omega;
	% index the orthogonalized sample matrix for relevant blocks
	U2 = orth(Y(1:n/2,1:k)); U3 = orth(Y((n/2)+1:end,k+1:end));
	% allocate space for a 'U' matrix
	U_mtrx = zeros(n,2*k);
	% interlace orthonormal U matrix
	U_mtrx(1:n/2,k+1:end) = U2; U_mtrx((n/2)+1:end,1:k) = U3;
	% create the Z matrix
	Z_mtrx = K_mtrx' * U_mtrx;
	% retrieve the relevant blocks 
	Z2 = Z_mtrx(1:n/2,1:k); Z3 = Z_mtrx((n/2)+1:end,k+1:end);
	% place in tree.
	% K2
	U_tree = U_tree.set(2,U2); Z_tree = Z_tree.set(2,Z2);
	% K3
	U_tree = U_tree.set(3,U3); Z_tree = Z_tree.set(3,Z3);
	end

	function output = ComputeNextLayer(l)
	% create a random matrix and interlace it with zeros.
	% idx_c: index cell array.
	idx_c = TraverseLayer(idx_tree,l);
	omega = randn(n,2*k); omega = InterlaceMtrx(omega,idx_c);
	% create sample matrix
	Y = K_mtrx * omega;
	% preallocate for upcoming result.
	U_mtrx = zeros(n,2*k);
%	MakeUMtrx(U_mtrx,Y,idx_c);
	for idx = 2:2:length(idx_c) 
		[s1 f1 s2 f2] = GetIntervals(idx_c,idx);
% TODO:
% define U1 and U2
% Find place to store -> the three.
% Define ComputeResidual
% Compute Z1 and Z2
		U1()
		B = ComputeResidual(); U1 = orth(U1 - B);
		B = ComputeResidual(); U2 = orth(U2 - B);
		U_mtrx(s2:f2,k+1:end) = U1; U_mtrx(s1:f1,1:k) = U2;
	end
%	MakeZMtrx(Z_mtrx,U,idx_c);
	for idx = 2:2:length(idx_c) 
		[s1 f1 s2 f2] = GetIntervals(idx_c,idx);
		B = ComputeResidual(); Z1 = U1 - B;
		B = ComputeResidual(); Z2 = U2 - B;
	end

	end

	% notice that this is in reverse chronological order.
	function residual = ComputeResidualU(m)
		for idx = 2:2:length(idx_c)
			B = table(U_tree.get(parent_node)).Var1(1:m/2)...
				* (table(Z_tree.get(parent_node+1)).Var1(m/2+1:end)'...
				* omega()); 
			if 
				B = B + ComputeResidual(m);
			end
			B = table(U_tree.get(parent_node)).Var1(m/2+1,end)...
				* (table(Z_tree.get(parent_node+1)).Var1(1:m/2)'...
				* omega()); 
			if
				B = B + ComputeResidual();
			end
		end
	end

	function interlaced_mtrx = InterlaceMtrx(mtrx,idx_c)
		for idx = 2:2:length(idx_c) 
			[s1 f1 s2 f2] = GetIntervals(idx_c,idx);
			omega(s1:f1,1:k) = 0; omega(s2:f2,k+1:end) = 0;
		end
	end

	function [U1 U2] = PartitionMtrx(mtrx,idx_c,idx)
		[s1 f1 s2 f2] = GetIntervals(idx_c,idx);
		U1 = Y(s1:f1,1:k); U2 = Y(s2:f2,k+1:end);
	end

	function [s1 s2 f1 s2] = GetIntervals(idx_c,idx)
		s1 = table(idx_c{idx-1}).Var1(1); f1 = table(idx_c{idx-1}).Var1(2);
		s2 = table(idx_c{idx}).Var1(1); f2 = table(idx_c{idx}).Var1(2);
	end

	function idx_c = TraverseLayer(btree,l)
	% ARGS:
	%	btree: binary tree with at least 1 layer.
	%	l: the layer whose nodes must be retrieved.
	% OUTPUT: cell array containing nxn array elements.
	% DESCRIPTION: loop through layer l's nodes and add them to cell array.
		it = btree.breadthfirstiterator; s = 2^l; f = 2*s-1;
		[idx_c{1:(f-s)}] = deal([]);
		for idx = s:f	idx_c{idx-(s-1)} = btree.get(it(idx));	end
	end
end
