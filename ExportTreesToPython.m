% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
n = 1000; diag_size = 125; k = 10; I = [1 n^d];
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
K_mtrx = rbf(ep,DM);
[U_tree Z_tree idx_tree] = NumericalExperiment(K_mtrx,n^d,k,diag_size,I);
exportToPython(U_tree,k,'U_tree')
exportToPython(Z_tree,k,'Z_tree')

function exportToPython(tree,k,name)
	itArr = tree.breadthfirstiterator;
	nodeLength = length(tree.get(itArr(2)));
	treeContainer = zeros(nodeLength,k,length(itArr));
	for iNode=2:length(itArr)
		nodeLength = length(tree.get(itArr(iNode)));
		treeContainer(1:nodeLength,:,itArr(iNode)) = tree.get(itArr(iNode));
	end
	save(sprintf('../HODLR_SLQ/%s',name),'treeContainer','itArr','-v7')
end
