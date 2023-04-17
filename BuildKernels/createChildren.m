function tree = createChildren(tree,node,diag_size)
	% subdivide at the midpoint; can be changed for more partitions
	interval = tree.get(node);
	interval_length = interval(2) - interval(1);
	midpoint = ceil(mean(interval));
	if interval_length < diag_size
		[tree leaf] = tree.addnode(node,interval);
		return
	end
	[tree left_node] = tree.addnode(node,[interval(1) midpoint-1]);
	[tree right_node] = tree.addnode(node,[midpoint interval(2)]);
	tree = createChildren(tree,left_node,diag_size);
	tree = createChildren(tree,right_node,diag_size);
end
