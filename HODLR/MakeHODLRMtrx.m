function output = MakeHODLRMtrx(kMtrxFcn,n,k,diagSize,I)
	% initialize index tree.
	idxTree = tree(I);
	idxTree = createChildren(idxTree,1,diagSize);
	uTree = idxTree; zTree = idxTree;
	% main
	computeFirstGeneration();
	treeDepth = floor(log2(nnodes(idxTree)+1));
	for g = 2:treeDepth-1
		computeNextGeneration(g);
	end
	leavesCell = computeLeaves(treeDepth);
	output = {uTree,zTree,leavesCell,idxTree};

	% helper functions.

	function computeFirstGeneration()
		% reference: PG. Martinsson, J. Tropp, 2020: pg 100-101
		% randomly sample the largest blocks of the matrix.
		omega = randn(n,2*k);
		omega(1:n/2,1:k) = 0; omega((n/2)+1:end,k+1:end) = 0;
		Y = kMtrxFcn(omega);
		U2 = orth(Y(1:n/2,1:k)); U3 = orth(Y((n/2)+1:end,k+1:end));
		uMtrx = zeros(n,2*k);
		% interlace in the opposite order (Fig 13).
		uMtrx(1:n/2,k+1:end) = U2; uMtrx((n/2)+1:end,1:k) = U3;
        % need an optional input that allows for adjoint
		zMtrx = kMtrxFcn(uMtrx);
		Z2 = zMtrx(1:n/2,1:k); Z3 = zMtrx((n/2)+1:end,k+1:end);
		uTree = uTree.set(2,U2); zTree = zTree.set(2,Z2);
		uTree = uTree.set(3,U3); zTree = zTree.set(3,Z3);
	end

	function computeNextGeneration(generation)
		% each generation has associated indices in the matrix; these indices
		% are organized in an index tree.
		% The zTree and uTree have associated cell arrays for the compressed
		% versions of the matrix blocks.
		[startIdx,idxCell] = traverseGeneration(idxTree,generation);
		itIdx = startIdx-1;
		it = idxTree.breadthfirstiterator;
		omega = randn(n,2*k);
		omega = interlaceMtrx(omega,idxCell);
		Y = kMtrxFcn(omega);
		uMtrx = zeros(n,2*k);
		for idx = 2:2:length(idxCell)
			[s1,f1,s2,f2] = getIntervals(idxCell,idx);
			U1 = Y(s1:f1,1:k);
			U2 = Y(s2:f2,k+1:end);
			U1 = orth(U1 - computeResidual(uTree,zTree,[s1 f1],omega,0));
			U2 = orth(U2 - computeResidual(uTree,zTree,[s2 f2],omega,0));
			uMtrx(s1:f1,k+1:end) = U1;
			uMtrx(s2:f2,1:k) = U2;
			itIdx = itIdx+2;
			uTree = uTree.set(it(itIdx-1),U1);
			uTree = uTree.set(it(itIdx),U2);
		end
		zMtrx = kMtrxFcn(uMtrx);
		itIdx = startIdx-1;
		for idx = 2:2:length(idxCell) 
			[s1,f1,s2,f2] = getIntervals(idxCell,idx);
			Z1 = zMtrx(s1:f1,1:k); Z2 = zMtrx(s2:f2,k+1:end);
			Z1 = Z1 - computeResidual(zTree,uTree,[s1 f1],uMtrx,0);
			Z2 = Z2 - computeResidual(zTree,uTree,[s2 f2],uMtrx,0);
			itIdx = itIdx+2;
			zTree = zTree.set(it(itIdx-1),Z1);
			zTree = zTree.set(it(itIdx),Z2);
		end
	end

	function leavesCell = computeLeaves(generation)
		% stack (number of leaves) identity blocks:
		idxCell = traverseLeafGeneration(idxTree,generation);
		eyeBlockLen = table(idxCell{1}).Var1(2);
		eyeBlocksNum = length(idxCell);
		% allocate space for leaves
		[leavesCell{1:eyeBlocksNum}] = deal([]);
		eyeStack = repmat(eye(eyeBlockLen),eyeBlocksNum,1);
		output = kMtrxFcn(eyeStack);
		for idx = 2:2:eyeBlocksNum
			[s1,f1,s2,f2] = getIntervals(idxCell,idx);
			leaf1 = output(s1:f1,:)-...
					computeResidual(uTree,zTree,[s1 f1],eyeStack,1);
			leaf2 = output(s2:f2,:)-...
					computeResidual(uTree,zTree,[s2 f2],eyeStack,1);
			leavesCell{idx-1} = leaf1;
			leavesCell{idx} = leaf2;
		end
	end

	function B = computeResidual(tree1,tree2,currInterval,randMtrx,isLeaf)
		nodePath = findNodeDFS(currInterval);
		startIdx = 2;
		if isLeaf
			startIdx = 1;
		end
		B = 0;
		for idx = startIdx:length(nodePath)
			[mtrx1,mtrx2,prevInterval] = getAncestor(tree1,tree2,...
													nodePath(idx));
			[s1,f1,s2,f2] = getAncestorIntervals(mtrx1,currInterval,...
												prevInterval);
			s3 = 1;
			f3 = table(size(randMtrx)).Var1(2);
			if ~isLeaf && isUpperTriangleBlock(nodePath(1))
				f3 = k;
			elseif ~isLeaf
				s3 = k+1;
			end
			B = B+mtrx1(s1:f1,:)*(mtrx2'*randMtrx(s2:f2,s3:f3));
		end
	end
% this function has bugs. f1 == 0 is a bad choice; should maybe be f1==1
	function [s1,f1,s2,f2] = getAncestorIntervals(mtrx,currInterval,prevInterval)
		s1 = mod(currInterval(1),length(mtrx));
		f1 = mod(currInterval(2),length(mtrx));
        %if s1 == 0, s1 = 1; end
		if f1 == 0, f1 = length(mtrx);	end
		% the columns of the previous block correspond to the rows of the
		% skinny matrix, by definition of matrix multiplication.
		s2 = prevInterval(1); f2 = prevInterval(2);
	end

	function [mtrx1,mtrx2,interval] = getAncestor(tree1,tree2,blockNum)
		if isUpperTriangleBlock(blockNum)
			newBlockNum = blockNum+1;
		else
			newBlockNum = blockNum-1;
		end
		mtrx1 = tree1.get(blockNum); mtrx2 = tree2.get(newBlockNum);
		interval = idxTree.get(newBlockNum);
	end

	function boolVal = isUpperTriangleBlock(blockNum)
		boolVal = (mod(blockNum,2) == 0);
	end

	function nodePath = findNodeDFS(interval)
		% find a node using depth-first search.
		itArr = idxTree.depthfirstiterator;
		for idx = 1:length(itArr)
			if all(idxTree.get(itArr(idx)) == interval)
				nodePath = idxTree.findpath(itArr(idx), 1);
				nodePath = nodePath(1:end-1);
				return
			end
		end
	end

	function mtrx = interlaceMtrx(mtrx,idxCell)
		for idx = 2:2:length(idxCell) 
			[s1,f1,s2,f2] = getIntervals(idxCell,idx);
			mtrx(s1:f1,1:k) = 0; mtrx(s2:f2,k+1:end) = 0;
		end
	end

	function [s1,f1,s2,f2] = getIntervals(idxCell,idx)
		s1 = table(idxCell{idx-1}).Var1(1);
		f1 = table(idxCell{idx-1}).Var1(2);
		s2 = table(idxCell{idx}).Var1(1);
		f2 = table(idxCell{idx}).Var1(2);
	end

	function [s,idxCell] = traverseGeneration(btree,g)
		it = btree.breadthfirstiterator;
		s = 2^g;
		f = 2*s-1;
		[idxCell{1:(f-s)}] = deal([]);
		for idx = s:f
			idxCell{idx-(s-1)} = btree.get(it(idx));
		end
	end

	function idxCell = traverseLeafGeneration(btree,g)
		it = btree.breadthfirstiterator;
		s = 2^g;
		f = s + s/2-1;
		[idxCell{1:(f-s)}] = deal([]);
		for idx=s:f
			idxCell{idx-(s-1)} = btree.get(it(idx));
		end
	end
end
