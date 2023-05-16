function y = HODLRMatVec(A,b)
% To test, try a transpose on U and Z.
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
	for ii = 2:2:nonLeafNodes
        % top block
		interval = idxTree.get(it(ii));
		s1 = interval(1);
		f1 = interval(2);
		u1 = uTree.get(it(ii));
		z1 = zTree.get(it(ii+1));
        % bottom block
        interval = idxTree.get(it(ii+1));     
		s2 = interval(1);
		f2 = interval(2);
		u2 = uTree.get(it(ii+1));
		z2 = zTree.get(it(ii));
		y(s1:f1,:) = y(s1:f1,:) + u1*(z1'*b(s2:f2,:));
        y(s2:f2,:) = y(s2:f2,:) + u2*(z2'*b(s1:f1,:));

	end
	% leaves 
    for ii = nonLeafNodes+1:numNodes
		interval = idxTree.get(it(ii));
		s = interval(1);
        f = interval(2);
		y(s:f,:) = y(s:f,:) + leavesCell{ii-nonLeafNodes}*b(s:f,:);
    end
end