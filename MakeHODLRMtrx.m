function output = MakeHODLRMtrx(Kmtrx,n,k,diagSize,I)
	% make binary tree
	idxTree = tree(I);
	idxTree = createChildren(idxTree,1,diagSize);
	% store U and Z values in nodes respective to index intervals.
	uTree = idxTree; zTree = idxTree;

	% Layer 1
	computeFirstGeneration();

	% we need to subtract away residuals for every generation after generation
	% 1.
	treeDepth = log2(nnodes(idxTree)+1) - 1;
	% subsequent generations. 
	for g = 2:treeDepth
		% compute the U and Z skinny matrices emerging from the next
		% generation's indices
		computeNextGeneration(g);
	end
	output = [uTree zTree idxTree];

	function output = computeFirstGeneration()
		% create a random matrix and interlace it with zeroes.
		omega = randn(n,2*k);
		omega(1:n/2,1:k) = 0; omega((n/2)+1:end,k+1:end) = 0;
		% create sample matrix
		Y = Kmtrx * omega;
		% index the orthogonalized sample matrix for relevant blocks
		U2 = orth(Y(1:n/2,1:k)); U3 = orth(Y((n/2)+1:end,k+1:end));
		% allocate space for a 'U' matrix
		uMtrx = zeros(n,2*k);
		% interlace orthonormal U matrix
		uMtrx(1:n/2,k+1:end) = U2; uMtrx((n/2)+1:end,1:k) = U3;
		% create the Z matrix
		zMtrx = Kmtrx' * uMtrx;
		% retrieve the relevant blocks 
		Z2 = zMtrx(1:n/2,1:k); Z3 = zMtrx((n/2)+1:end,k+1:end);
		% place in tree.
		% K2,3 -> K2; enumerate based on the first index.
		uTree = uTree.set(2,U2); zTree = zTree.set(2,Z2);
		% K3
		uTree = uTree.set(3,U3); zTree = zTree.set(3,Z3);
		end

	function noReturn = computeNextGeneration(generation)
		% ARGS:
		%	generation: integer representing the relevant generation in the
		%	tree.
		%	it: iterator representing to first node at this generation in the
		%	tree.
		% OUTPUT: no output; assign values to this generation's U and Z trees.
		% DESCRIPTION: perform the peeling algorithm steps described in
		% Martinsson/Tropp Section 20.5 (pgs 98-101).

		% idxCell: index cell array of this generation's index pairs; all nodes
		% in single generation of 'idxTree'
		[sIdx idxCell] = traverseGeneration(idxTree,generation);
		itIdx = sIdx-1; it = idxTree.breadthfirstiterator;
		% create a random matrix and interlace it with zeros.
		omega = randn(n,2*k); omega = interlaceMtrx(omega,idxCell);
		% create sample matrix Y.
		Y = Kmtrx * omega;
		% preallocate for upcoming result.
		uMtrx = zeros(n,2*k);
		% Create U vectors and interlace them with zeros.
		for idx = 2:2:length(idxCell)
			% determine which blocks in the sample matrix we wish to index.
			[s1 f1 s2 f2] = getIntervals(idxCell,idx);
			U1 = Y(s1:f1,1:k); U2 = Y(s2:f2,k+1:end);

			% Compute residuals and orthogonalize the 'purified' result.
			U1 = orth(U1 - computeResidualU([s1 f1],omega));
			U2 = orth(U2 - computeResidualU([s2 f2],omega));
			% interlace U1 and U2 into relevant blocks.
			uMtrx(s1:f1,k+1:end) = U1; uMtrx(s2:f2,1:k) = U2;
			% set U values and adjust the iterator.
			itIdx = itIdx+2;
			uTree = uTree.set(it(itIdx),U2);
			uTree = uTree.set(it(itIdx-1),U1);
		end
		% create sample Z matrix
		zMtrx = Kmtrx' * uMtrx;
		% Create U vectors and interlace them with zeros.
		itIdx = sIdx-1;
		for idx = 2:2:length(idxCell) 
			% determine which blocks in the sample matrix we wish to index.
			[s1 f1 s2 f2] = getIntervals(idxCell,idx);
			Z1 = zMtrx(s1:f1,1:k); Z2 = zMtrx(s2:f2,k+1:end);
			% Compute residuals and orthogonalize the 'purified' result.
			Z1 = Z1 - computeResidualZ([s1 f1],uMtrx);
			Z2 = Z2 - computeResidualZ([s2 f2],uMtrx);
			% set Z values and adjust the iterator.
			itIdx = itIdx+2;
			zTree = zTree.set(it(itIdx),Z2);
			zTree = zTree.set(it(itIdx-1),Z1);
		end
	end

	function B = computeResidualU(currInterval,omega)
		% ARGS:
		%	currInterval: array containing start and finish for sample matrix
		%	block.
		%   omega: random skinny matrix interlaced with zeros.
		% OUTPUT: vector containing residual value for sample matrix block.
		% DESCRIPTION: (1) find the array matching 'currInterval' and (2) the
		% path of parent nodes. Compute the residual by (3) matching every
		% parent with a portion of the sample matrix.

		% (1) and (2)
		nodePath = findNodeDFS(currInterval);
		B = 0;
		% (3)
		% start with 2nd index; first index is node itself.
		for idx = 2:length(nodePath)
			% retrieve 'idx'th ancestor.
			% prevInterval comes from the sibling of the ancestor node. It
			% determines which rows of the sample matrix to index.
			[U Z prevInterval] = getAncestorU(idx,nodePath);
			% what portion of the 'ancestor block' to index is determined by
			% which rows the the current interval is in.
			s1 = mod(currInterval(1),length(U));
			f1 = mod(currInterval(2),length(U));
			if f1 == 0	f1 = length(U);	end

			% the columns of the previous block correspond to the rows of the
			% skinny matrix, by definition of matrix multiplication.
			s2 = prevInterval(1); f2 = prevInterval(2);

			% if a block is in the upper triangle, then we index the first k
			% columns of omega.

			if isUpperTriangleBlock(nodePath(1))
				B = B + U(s1:f1,:) * (Z' * omega(s2:f2,1:k));
			else
				B = B + U(s1:f1,:) * (Z' * omega(s2:f2,k+1:end));
			end
		end
	end

	function B = computeResidualZ(currInterval,uMtrx)
		% ARGS:
		%	interval: array containing start and finish for sample matrix
		%	block.
		%	uMtrx: zero matrix interlaced with U skinny sub-matrices.
		% OUTPUT: vector containing residual value for sample matrix block.
		% DESCRIPTION: (1) find the array matching 'currInterval' and (2) the
		% path of parent nodes. Compute the residual by (3) matching every
		% parent with a portion of the sample matrix.

		% (1) and (2)
		nodePath = findNodeDFS(currInterval);
		B = 0;
		% (3)
		for idx = 2:length(nodePath)
			% retrieve 'idx'th ancestor.
			[U Z prevInterval] = getAncestorZ(idx,nodePath);
			% what portion of the 'ancestor block' to index is determined by
			% which rows the the current interval is in.
			s1 = mod(currInterval(1),length(Z));
			f1 = mod(currInterval(2),length(Z));
			if f1 == 0	f1 = length(Z);	end

			% the columns of the previous block correspond to the rows of the
			% skinny matrix, by definition of matrix multiplication.
			s2 = prevInterval(1); f2 = prevInterval(2);

			% if a block is in the upper triangle, then we index the first k
			% columns of uMtrx.
			if isUpperTriangleBlock(nodePath(1))
				B = B + Z(s1:f1,:) * (U' * uMtrx(s2:f2,1:k));
			else
				B = B + Z(s1:f1,:) * (U' * uMtrx(s2:f2,k+1:end));
			end
		end
	end

	% Node checks out.
	function [U Z interval] = getAncestorU(idx,nodePath)
		% ARGS:
		%	idx: index of the nodePath indicating which ancestor to retrieve.
		%	nodePath: path of ancestor nodes in the idxTree.
		% OUTPUT: array containing an ancestor node's corresponding U block, Z
		% block, and interval.
		% DESCRIPTION: depending on whether the node corresponds to a block in
		% the upper or lower triangle, the sibling is either to the left or
		% right of the index in question.

		% take the first index of the block in the Tropp diagram. That
		% corresponds to the node.
		blockNum = nodePath(idx);
		if isUpperTriangleBlock(blockNum)
			U = uTree.get(blockNum); Z = zTree.get(blockNum+1);
			% the interval of the Z node is the interval of the columns of
			% sample matrix.
			interval = idxTree.get(blockNum+1);
		else
			U = uTree.get(blockNum); Z = zTree.get(blockNum-1);
			% the interval of the Z node is the interval of the columns of
			% sample matrix.
			interval = idxTree.get(blockNum-1);
		end
	end

	function [U Z interval] = getAncestorZ(idx,nodePath)
		% ARGS:
		%	idx: index of the nodePath indicating which ancestor to retrieve.
		%	nodePath: path of ancestor nodes in the idxTree.
		% OUTPUT: array containing an ancestor node's corresponding U block, Z
		% block, and interval.
		% DESCRIPTION: Same process as the for U, but the order of which node
		% is pulled from which tree is reversed.

		blockNum = nodePath(idx);
		if isUpperTriangleBlock(blockNum)
			% order is reversed.
			Z = zTree.get(blockNum); U = uTree.get(blockNum+1);
			% the interval of the Z node is the interval of the columns of
			% sample matrix.
			interval = idxTree.get(blockNum+1);
		else
			% order is reversed.
			Z = zTree.get(blockNum); U = uTree.get(blockNum-1);
			% the interval of the Z node is the interval of the columns of
			% sample matrix.
			interval = idxTree.get(blockNum-1);
		end
	end

	function boolVal = isUpperTriangleBlock(blockNum)
		% ARGS:
		%	blockNum: integer value corresponding to submatrix blocks/node
		%	number.
		% OUTPUT: boolean value.
		% DESCRIPTION: all blocks in the upper triangle are defined to have
		% even values.
		boolVal = (mod(blockNum,2) == 0);
	end

	function nodePath = findNodeDFS(interval)
		% ARGS:
		%	interval: array with indices corresponding to a node in idxTree.
		% OUTPUT: array of all ancestors of a node containing 'interval'.
		% DESCRIPTION: do a depth first search to retrieve all of a given
		% node's ancestors.
		
		% iterator array with DFS arrangement of indices.
		itArr = idxTree.depthfirstiterator;
		for idx = 1:length(itArr)
			% find the node in the tree by comparing to other intervals in
			% idxTree.
			if all(idxTree.get(itArr(idx)) == interval)
				% get all of the ancestor nodes leading back to index 1.
				nodePath = idxTree.findpath(itArr(idx), 1);
				% discard the first (irrelevant) and last node (node can't be a
				% parent of itself).
				nodePath = nodePath(1:end-1);
				return
			end
		end
	end

	function mtrx = interlaceMtrx(mtrx,idxCell)
		% ARGS:
		%	mtrx: input matrix to be interlaced.
		%	idxCell: cell array containing indices at this genration in the tree.
		% OUTPUT: same input matrix interlaced with zeros in some blocks.
		% DESCRIPTION: loop through each block corresponding to the intervals
		% at a given generation and set that block to zero.
		for idx = 2:2:length(idxCell) 
			[s1 f1 s2 f2] = getIntervals(idxCell,idx);
			mtrx(s1:f1,1:k) = 0; mtrx(s2:f2,k+1:end) = 0;
		end
	end

	function [s1 f1 s2 f2] = getIntervals(idxCell,idx)
		% ARGS:
		% 	idxCell: cell array with nodes of the idxTree for this generation.
		%	idx: node in the idxTree whose start and finish values we need.
		% OUTPUT: array containing the start and finish indices of two sibling
		% nodes in idxTree.
		% DESCRIPTION: Extract the start and finish indices of two sibling
		% nodes in the idxTree by traversing those nodes through idxCell.

		% start and finish indices of sibling
		s1 = table(idxCell{idx-1}).Var1(1); f1 = table(idxCell{idx-1}).Var1(2);
		% start and finish indices of node represented by 'idx'.
		s2 = table(idxCell{idx}).Var1(1); f2 = table(idxCell{idx}).Var1(2);
	end

	function [s idxCell] = traverseGeneration(btree,g)
		% ARGS:
		%	btree: binary tree with at least 1 generation.
		%	g: the generation whose nodes must be retrieved.
		% OUTPUT: cell array containing a 2x1 array in each cell index.
		% DESCRIPTION: loop through generation nodes and add them to cell
		% array.
		it = btree.breadthfirstiterator; s = 2^g; f = 2*s-1;
		% preallocate memory.
		[idxCell{1:(f-s)}] = deal([]);
		for idx = s:f	idxCell{idx-(s-1)} = btree.get(it(idx));	end
	end
end
