function y = HODLRMatVec(A,b)
	uTree = A{1};
	zTree = A{2};
	leavesCell = A{3};
	idxTree = A{4};

	% preallocate y; same size as b since A is square.
	y = zeros(size(b));

	% determine the number of nodes.
	numNodes = nnodes(idxTree);
	% determine number of non-leaf nodes
	nonLeafNodes = numNodes-length(leavesCell);

	% set up iterator.
	it = idxTree.breadthfirstiterator;

	% traverse all nodes of the binary trees in pairs, except the last leaves
	% layer.
	for ii = 2:nonLeafNodes
		interval = idxTree.get(it(ii));
		s = interval(1); f = interval(2);
		u = uTree.get(it(ii));
		z = zTree.get(it(ii+1));
		if mod(ii,2) ~= 0
			u = uTree.get(it(ii));
			z = zTree.get(it(ii-1));
		end
		y(s:f) = y(s:f) + u*(z'*b(s:f));
	end
	% leaves
	for ii = nonLeafNodes+1:numNodes
		interval = idxTree.get(it(ii));
		s = interval(1); f = interval(2);
		y(s:f) = y(s:f) + leavesCell{ii-nonLeafNodes}*b(s:f);
	end
end