function [U_tree Z_tree idx_tree] = NumericalExperiment(K_mtrx,n,k,diag_size,I)
	% make binary tree
	idx_tree = tree(I);
	idx_tree = createChildren(idx_tree,1,diag_size);
	% store U and Z values in nodes respective to index intervals.
	U_tree = idx_tree; Z_tree = idx_tree;

	% Layer 1
	computeFirstGeneration();

	% we need to subtract away residuals for every generation after generation
	% 1.
	tree_depth = log2(nnodes(idx_tree)+1) - 1;
	% subsequent generations. 
	for g = 2:tree_depth
		% compute the U and Z skinny matrices emerging from the next
		% generation's indices
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
		% K2,3 -> K2; enumerate based on the first index.
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

		% idx_c: index cell array of this generation's index pairs; all nodes
		% in single generation of 'idx_tree'
		[starting_idx idx_c] = traverseGeneration(idx_tree,generation);
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
			U1 = Y(s1:f1,1:k); U2 = Y(s2:f2,k+1:end);

			% Compute residuals and orthogonalize the 'purified' result.
			U1 = orth(U1 - computeResidualU([s1 f1],omega));
			U2 = orth(U2 - computeResidualU([s2 f2],omega));
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
			Z1 = Z_mtrx(s1:f1,1:k); Z2 = Z_mtrx(s2:f2,k+1:end);
			% Compute residuals and orthogonalize the 'purified' result.
			Z1 = Z1 - computeResidualZ([s1 f1],U_mtrx);
			Z2 = Z2 - computeResidualZ([s2 f2],U_mtrx);
			% set Z values and adjust the iterator.
			it_idx = it_idx+2;
			Z_tree = Z_tree.set(it(it_idx),Z2);
			Z_tree = Z_tree.set(it(it_idx-1),Z1);
		end
	end

	function B = computeResidualU(curr_interval,omega)
		% ARGS:
		%	curr_interval: array containing start and finish for sample matrix
		%	block.
		%   omega: random skinny matrix interlaced with zeros.
		% OUTPUT: vector containing residual value for sample matrix block.
		% DESCRIPTION: (1) find the array matching 'curr_interval' and (2) the
		% path of parent nodes. Compute the residual by (3) matching every
		% parent with a portion of the sample matrix.

		% (1) and (2)
		node_path = findNodeDFS(curr_interval);
		B = 0;
		% (3)
		% start with 2nd index; first index is node itself.
		for idx = 2:length(node_path)
			% retrieve 'idx'th ancestor.
			% prev_interval comes from the sibling of the ancestor node. It
			% determines which rows of the sample matrix to index.
			[U Z prev_interval] = getAncestorU(idx,node_path);
			% what portion of the 'ancestor block' to index is determined by
			% which rows the the current interval is in.
			s1 = mod(curr_interval(1),length(U));
			f1 = mod(curr_interval(2),length(U));
			if f1 == 0	f1 = length(U);	end

			% the columns of the previous block correspond to the rows of the
			% skinny matrix, by definition of matrix multiplication.
			s2 = prev_interval(1); f2 = prev_interval(2);

			% if a block is in the upper triangle, then we index the first k
			% columns of omega.

			if isUpperTriangleBlock(node_path(1))
				B = B + U(s1:f1,:) * (Z' * omega(s2:f2,1:k));
			else
				B = B + U(s1:f1,:) * (Z' * omega(s2:f2,k+1:end));
			end
		end
	end

	function B = computeResidualZ(curr_interval,U_mtrx)
		% ARGS:
		%	interval: array containing start and finish for sample matrix
		%	block.
		%	U_mtrx: zero matrix interlaced with U skinny sub-matrices.
		% OUTPUT: vector containing residual value for sample matrix block.
		% DESCRIPTION: (1) find the array matching 'curr_interval' and (2) the
		% path of parent nodes. Compute the residual by (3) matching every
		% parent with a portion of the sample matrix.

		% (1) and (2)
		node_path = findNodeDFS(curr_interval);
		B = 0;
		% (3)
		for idx = 2:length(node_path)
			% retrieve 'idx'th ancestor.
			[U Z prev_interval] = getAncestorZ(idx,node_path);
			% what portion of the 'ancestor block' to index is determined by
			% which rows the the current interval is in.
			s1 = mod(curr_interval(1),length(Z));
			f1 = mod(curr_interval(2),length(Z));
			if f1 == 0	f1 = length(Z);	end

			% the columns of the previous block correspond to the rows of the
			% skinny matrix, by definition of matrix multiplication.
			s2 = prev_interval(1); f2 = prev_interval(2);

			% if a block is in the upper triangle, then we index the first k
			% columns of U_mtrx.
			if isUpperTriangleBlock(node_path(1))
				B = B + Z(s1:f1,:) * (U' * U_mtrx(s2:f2,1:k));
			else
				B = B + Z(s1:f1,:) * (U' * U_mtrx(s2:f2,k+1:end));
			end
		end
	end

	% Node checks out.
	function [U Z interval] = getAncestorU(idx,node_path)
		% ARGS:
		%	idx: index of the node_path indicating which ancestor to retrieve.
		%	node_path: path of ancestor nodes in the idx_tree.
		% OUTPUT: array containing an ancestor node's corresponding U block, Z
		% block, and interval.
		% DESCRIPTION: depending on whether the node corresponds to a block in
		% the upper or lower triangle, the sibling is either to the left or
		% right of the index in question.

		% take the first index of the block in the Tropp diagram. That
		% corresponds to the node.
		block_num = node_path(idx);
		if isUpperTriangleBlock(block_num)
			U = U_tree.get(block_num); Z = Z_tree.get(block_num+1);
			% the interval of the Z node is the interval of the columns of
			% sample matrix.
			interval = idx_tree.get(block_num+1);
		else
			U = U_tree.get(block_num); Z = Z_tree.get(block_num-1);
			% the interval of the Z node is the interval of the columns of
			% sample matrix.
			interval = idx_tree.get(block_num-1);
		end
	end

	function [U Z interval] = getAncestorZ(idx,node_path)
		% ARGS:
		%	idx: index of the node_path indicating which ancestor to retrieve.
		%	node_path: path of ancestor nodes in the idx_tree.
		% OUTPUT: array containing an ancestor node's corresponding U block, Z
		% block, and interval.
		% DESCRIPTION: Same process as the for U, but the order of which node
		% is pulled from which tree is reversed.

		block_num = node_path(idx);
		if isUpperTriangleBlock(block_num)
			% order is reversed.
			Z = Z_tree.get(block_num); U = U_tree.get(block_num+1);
			% the interval of the Z node is the interval of the columns of
			% sample matrix.
			interval = idx_tree.get(block_num+1);
		else
			% order is reversed.
			Z = Z_tree.get(block_num); U = U_tree.get(block_num-1);
			% the interval of the Z node is the interval of the columns of
			% sample matrix.
			interval = idx_tree.get(block_num-1);
		end
	end

	function bool_value = isUpperTriangleBlock(block_num)
		% ARGS:
		%	block_num: integer value corresponding to submatrix blocks/node
		%	number.
		% OUTPUT: boolean value.
		% DESCRIPTION: all blocks in the upper triangle are defined to have
		% even values.
		bool_value = (mod(block_num,2) == 0);
	end

	function node_path = findNodeDFS(interval)
		% ARGS:
		%	interval: array with indices corresponding to a node in idx_tree.
		% OUTPUT: array of all ancestors of a node containing 'interval'.
		% DESCRIPTION: do a depth first search to retrieve all of a given
		% node's ancestors.
		
		% iterator array with DFS arrangement of indices.
		it_a = idx_tree.depthfirstiterator;
		for idx = 1:length(it_a)
			% find the node in the tree by comparing to other intervals in
			% idx_tree.
			if all(idx_tree.get(it_a(idx)) == interval)
				% get all of the ancestor nodes leading back to index 1.
				node_path = idx_tree.findpath(it_a(idx), 1);
				% discard the first (irrelevant) and last node (node can't be a
				% parent of itself).
				node_path = node_path(1:end-1);
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

	function [s1 f1 s2 f2] = getIntervals(idx_c,idx)
		% ARGS:
		% 	idx_c: cell array with nodes of the idx_tree for this generation.
		%	idx: node in the idx_tree whose start and finish values we need.
		% OUTPUT: array containing the start and finish indices of two sibling
		% nodes in idx_tree.
		% DESCRIPTION: Extract the start and finish indices of two sibling
		% nodes in the idx_tree by traversing those nodes through idx_c.

		% start and finish indices of sibling
		s1 = table(idx_c{idx-1}).Var1(1); f1 = table(idx_c{idx-1}).Var1(2);
		% start and finish indices of node represented by 'idx'.
		s2 = table(idx_c{idx}).Var1(1); f2 = table(idx_c{idx}).Var1(2);
	end

	function [s idx_c] = traverseGeneration(btree,g)
		% ARGS:
		%	btree: binary tree with at least 1 generation.
		%	g: the generation whose nodes must be retrieved.
		% OUTPUT: cell array containing a 2x1 array in each cell index.
		% DESCRIPTION: loop through generation nodes and add them to cell
		% array.
		it = btree.breadthfirstiterator; s = 2^g; f = 2*s-1;
		% preallocate memory.
		[idx_c{1:(f-s)}] = deal([]);
		for idx = s:f	idx_c{idx-(s-1)} = btree.get(it(idx));	end
	end
end
