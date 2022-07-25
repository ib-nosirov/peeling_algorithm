% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 20; d = 1;
n = 1000; diagSize= 125; k = 10; I = [1 n^d];
% Point evals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
Kmtrx = rbf(ep,DM);
b = ones(n,1);
%y = Kmtrx * b;
%y;

exportToPython(MakeHODLRMtrx(kMtrx,n^d,k,diagSize,I),'HODLR_mtrx.mat');

function noReturn = exportToPython(HODLRMtrx,name)
	u_tree = prepTree(HODLRMtrx{1});
	z_tree = prepTree(HODLRMtrx{2});
	leaves_cell = HODLRMtrx{3};
	idx_tree = prepTree(HODLRMtrx{4});
	save(sprintf('../HODLR_SLQ/%s',name),'u_tree','z_tree','leaves_cell','idx_tree','-v7');
end

function treeContainer = prepTree(tree)
	itArr = tree.breadthfirstiterator;
	% preallocate memory.
	% note that we are discarding the first index 
	[treeContainer{1:length(itArr)-1}] = deal([]);
	for iNode=2:length(itArr)
		treeContainer{iNode-1} = tree.get(itArr(iNode));
	end
end
