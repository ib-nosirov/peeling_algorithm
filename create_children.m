function tree = create_children(tree, node, diag_size)
    % subdivide at the midpoint; can be changed for more partitions
    interval = tree.get(node);
    interval_length = interval(2) - interval(1);
    if interval_length > diag_size
        midpoint = ceil(mean(interval));
	[tree left_node] = tree.addnode(node, [interval(1) ...
                                               midpoint-1]);
        [tree right_node] = tree.addnode(node, [midpoint ... 
                                                interval(2)]);
        tree = create_children(tree, left_node, diag_size);
        tree = create_children(tree, right_node, diag_size);
    end
end
