function output = HutchPlusPlus(K,U_tree,idx_tree)
	% ARGS:
	%	U_tree: k-d tree containing low-rank apporximations of submatrix A.
	%	idx_tree: k-d tree containing indices for the submatrices.
	output = 0;
	[N k] = size(U_tree.get(2)); tree_depth = log2(nnodes(idx_tree)+1) - 1;
	N = N * 2;
	for g = 1:tree_depth
		% Rademacher random variable.
		G = randi([0 1],N,2*k);

		% Get all of the U matrices from this layer

		% start by isolating the necessary indices.
		[starting_idx idx_c] = traverseGeneration(idx_tree,g);
		% preallocate for the large Q matrix.
		Q = zeros(N,2*k);
		it_idx = starting_idx; it = idx_tree.breadthfirstiterator;
		%  
		for idx = 2:2:length(idx_c)
			% determine which blocks in the sample matrix we wish to index.
			[s1 f1 s2 f2] = getIntervals(idx_c,idx);
			% interlace U1 and U2 into relevant blocks.
			Q(s1:f1,1:k) = U_tree.get(it(it_idx));
			Q(s2:f2,k+1:end) = U_tree.get(it(it_idx+1));
			it_idx = it_idx+2;
		end
		% Part 1: Compute tr(Q'*K*Q) exactly
		tracePart1 = trace(Q'*K*Q) % THIS USES 2r MULTIPLIES BY K
		tracePart2 = (1/size(G,2))*trace(G' * (eye(N)-Q*Q')*K*(eye(N)-Q*Q') * G)
		return
		output = output + tracePart1 + tracePart2;
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
end

