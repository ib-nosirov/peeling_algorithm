function [U_tree Z_tree idx_tree] = NumericalExperiment(K_mtrx,n,k,diag_size,I)
% make binary tree
idx_tree = tree(I);
idx_tree = createChildren(idx_tree,1,diag_size);
% store U and Z values in nodes respective to index intervals.
U_tree = idx_tree; Z_tree = idx_tree;

% Layer 1
computeFirstGeneration();

% we need to subtract away residuals for every generation after generation 1.
tree_depth = log2(nnodes(idx_tree)+1) - 1;
% subsequent generations. 
for g = 2:tree_depth
	computeNextGeneration(g);
end

	function output = computeFirstGeneration()
		% create a random matrix and interlace it with zeroes.
		omega = randn(n,2*k);
		omega(1:n/2,1:k) = 0; omega((n/2)+1:end,k+1:end) = 0;
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

	function no_return = computeNextGeneration(generation)
		% ARGS:
		%	generation: integer representing the relevant generation in the
		%	tree.
		%	it: iterator representing to first node at this generation in the
		%	tree.
		% OUTPUT: no output; assign values to this generation's U and Z trees.
		% DESCRIPTION: perform the peeling algorithm steps described in
		% Martinsson/Tropp Section 20.5 (pgs 98-101).

		% idx_c: index cell array.
		[starting_idx idx_c] = traverseLayer(idx_tree,generation);
		it_idx = starting_idx-1; it = idx_tree.breadthfirstiterator;
		% create a random matrix and interlace it with zeros.
		omega = randn(n,2*k); omega = interlaceMtrx(omega,idx_c);
		% create sample matrix Y.
		Y = K_mtrx * omega;
		% preallocate for upcoming result.
		U_mtrx = zeros(n,2*k);
		% Create U vectors and interlace them with zeros.
		for idx = 2:2:length(idx_c) 
			% determine which blocks in the sample matrix we wish to index.
			[s1 f1 s2 f2] = getIntervals(idx_c,idx);
			% Compute residuals and orthogonalize the 'purified' result.
			U1 = Y(s1:f1,1:k); U1 = orth(U1 - computeResidualU([s1 f1],omega));
			U2 = Y(s2:f2,k+1:end); U2 = orth(U2 - computeResidualU([s2 f2],...
										omega));
			% interlace U1 and U2 into relevant blocks.
			U_mtrx(s1:f1,k+1:end) = U1; U_mtrx(s2:f2,1:k) = U2;
			% set U values and adjust the iterator.
			it_idx = it_idx+2;
			U_tree = U_tree.set(it(it_idx),U2);
			U_tree = U_tree.set(it(it_idx-1),U1);
		end
		% create sample Z matrix
		Z_mtrx = K_mtrx' * U_mtrx;
		% Create U vectors and interlace them with zeros.
		it_idx = starting_idx-1;
		for idx = 2:2:length(idx_c) 
			% determine which blocks in the sample matrix we wish to index.
			[s1 f1 s2 f2] = getIntervals(idx_c,idx);
			% Compute residuals and orthogonalize the 'purified' result.
			Z1 = Z_mtrx(s1:f1,1:k); Z1 = Z1 - computeResidualZ([s1 f1],U_mtrx);
			Z2 = Z_mtrx(s2:f2,k+1:end); Z2 = Z2 - computeResidualZ([s2 f2],...
										U_mtrx);
			it_idx = it_idx+2;
			Z_tree = Z_tree.set(it(it_idx),Z2);
			Z_tree = Z_tree.set(it(it_idx-1),Z1);
		end

	end

	function B = computeResidualU(interval,omega)
		% ARGS:
		%	interval: array containing start and finish for sample matrix
		%	block.
		%   omega: random skinny matrix interlaced with zeros.
		% OUTPUT: vector containing residual value for sample matrix block.
		% DESCRIPTION: find the array matching 'interval', as well as the path
		% of parent nodes. Compute the residual by matching every parent with a
		% portion of the sample matrix.
		node_path = findNodeDFS(interval);
		B = 0;
		for idx = 1:length(node_path)
			[U Z node] = getPrevGeneration(idx,node_path);
			s1 = mod(interval(1),length(U));
			f1 = mod(interval(2),length(U));
			if f1 == 0	f1 = length(U);	end
			U = U(s1:f1,:);
			s2 = node(1); f2 = node(2);
			B = B + U * (Z' * omega(s2:f2,1:k));
		end
	end

	function B = computeResidualZ(interval,U_mtrx)
		% ARGS:
		%	interval: array containing start and finish for sample matrix
		%	block.
		%	U_mtrx: zero matrix interlaced with U skinny sub-matrices.
		% OUTPUT: vector containing residual value for sample matrix block.
		% DESCRIPTION: find the array matching 'interval', as well as the path
		% of parent nodes. Compute the residual by matching every parent with a
		% portion of the sample matrix.
		node_path = findNodeDFS(interval);
		B = 0;
		for idx = 1:length(node_path)
			[U Z node] = getPrevGeneration(idx,node_path);
			s1 = mod(interval(1),length(Z));
			f1 = mod(interval(2),length(Z));
			if f1 == 0	f1 = length(Z);	end
			Z = Z(s1:f1,:);
			s2 = node(1); f2 = node(2);
			B = B + Z * (U' * U_mtrx(s2:f2,1:k));
		end
	end

	% Node checks out.
	function [U Z node] = getPrevGeneration(idx,node_path)
		if isTopNode(node_path(idx))
			U = U_tree.get(node_path(idx));
			Z = Z_tree.get(node_path(idx)+1);
			node = idx_tree.get(node_path(idx)+1);
		else
			U = U_tree.get(node_path(idx)-1);
			Z = Z_tree.get(node_path(idx));
			node = idx_tree.get(node_path(idx)-1);
		end
	end

	% Helper method checks out.
	function truth_value = isTopNode(node_value)
		truth_value = (mod(node_value,2) == 0);
	end

	function node_path = findNodeDFS(interval)
		% iterator array with DFS arrangement of indices.
		it_a = idx_tree.depthfirstiterator;
		for idx = 1:length(it_a)
			if all(idx_tree.get(it_a(idx)) == interval)
				node_path = idx_tree.findpath(it_a(idx), 1);
				node_path = node_path(2:end-1);
				return
			end
		end
	end

	function mtrx = interlaceMtrx(mtrx,idx_c)
		% ARGS:
		%	mtrx: input matrix to be interlaced.
		%	idx_c: cell array containing indices at this genration in the tree.
		% OUTPUT: same input matrix interlaced with zeros in some blocks.
		% DESCRIPTION: loop through each block corresponding to the intervals
		% at a given generation and set that block to zero.
		for idx = 2:2:length(idx_c) 
			[s1 f1 s2 f2] = getIntervals(idx_c,idx);
			mtrx(s1:f1,1:k) = 0; mtrx(s2:f2,k+1:end) = 0;
		end
	end

	function [U1 U2] = partitionMtrx(mtrx,idx_c,idx)
		[s1 f1 s2 f2] = getIntervals(idx_c,idx);
		U1 = mtrx(s1:f1,1:k); U2 = mtrx(s2:f2,k+1:end);
	end

	function [s1 f1 s2 f2] = getIntervals(idx_c,idx)
		s1 = table(idx_c{idx-1}).Var1(1); f1 = table(idx_c{idx-1}).Var1(2);
		s2 = table(idx_c{idx}).Var1(1); f2 = table(idx_c{idx}).Var1(2);
	end

	function [s idx_c] = traverseLayer(btree,g)
		% ARGS:
		%	btree: binary tree with at least 1 generation.
		%	g: the generation whose nodes must be retrieved.
		% OUTPUT: cell array containing nxn array elements.
		% DESCRIPTION: loop through generation nodes and add them to cell
		% array.
		it = btree.breadthfirstiterator; s = 2^g; f = 2*s-1;
		[idx_c{1:(f-s)}] = deal([]);
		for idx = s:f	idx_c{idx-(s-1)} = btree.get(it(idx));	end
	end
end
