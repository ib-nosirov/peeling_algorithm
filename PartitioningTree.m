% Pseudocode: Create a binary tree that partitions until the
% difference between the elements of each node is smaller than the
% desired value n (typically in the few hundreds)

I = [0 1000];
root = create_node(I);
make_partition(root, 130);
output_tree(root);

% create node
% populate the node with a value
function node = create_node(interval)
    node = struct('left', [], 'right', [], 'value', interval);
end

% input: struct, output: [struct struct]
function [left_child right_child] = create_children(node)
    % subdivide at the midpoint; can be changed for more partitions
    midpoint = mean(node.value);
    left_child = create_node([node.value(1) midpoint]);
    right_child = create_node([midpoint node.value(2)]);
end

% recursively keep dividing until some threshold N has been reached.
function partition_tree = make_partition(node, diag_size)
    % base case: do nothing
    interval_length = node.value(2) - node.value(1);
    if diag_size < interval_length
         % new layer of binary tree
         [node.left node.right] = create_children(node);
         make_partition(node.left, diag_size);
         make_partition(node.right, diag_size);
    end 
end

% output tree
function node = output_tree(root)
    root
    output_tree(root.left);
    output_tree(root.right);
end
