% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
n = 1000; diagSize= 125; k = 10; I = [1 n^d];
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
Kmtrx = rbf(ep,DM);
exportToPython(MakeHODLRMtrx(Kmtrx,n^d,k,diagSize,I),k,'HODLR_mtrx.mat')

function noReturn = exportToPython(HODLRMtrx,approxLowRank,name)
	u_tree = prepTree(HODLRMtrx(1),approxLowRank);
	z_tree = prepTree(HODLRMtrx(2),approxLowRank);
	idx_tree = prepTree(HODLRMtrx(3),1);
	save(sprintf('../HODLR_SLQ/%s',name),'u_tree','z_tree','idx_tree','-v7');
end

function treeContainer = prepTree(tree,k)
	itArr = tree.breadthfirstiterator;
	nodeLength = length(tree.get(itArr(2)));
	% note that we are discarding the first index 
	treeContainer = zeros(nodeLength,k,length(itArr)-1);
	for iNode=2:length(itArr)
		nodeLength = length(tree.get(itArr(iNode)));
		treeContainer(1:nodeLength,:,iNode-1) = tree.get(itArr(iNode));
	end
end
