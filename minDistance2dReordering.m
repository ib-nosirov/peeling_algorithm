% Author: Ibrohim Nosirov
% Date: 2023-01-15
% Version: 1.0
% Description: a maximin algorithm for ordering points in 2 dimensions into a
% HODLR-type matrix. The algorithm is based on the idea that faraway
% interactions have near-zero, or even slightly negative values, and can
% therefore be compressed into HODLR-type if these interactions reside on the
% off-diagonal blocks of the matrix.
% Input: a set of points in 2 dimensions
% Output: a HODLR-type matrix

% Generate set of 2D points as an array of size 2xN.
% The points are generated randomly in the unit square.
clear;
close all;
N = 100;
% a set of (x,y) points with an appended 0/1 vector indicating whether each
% point has been visited.
dim2Points = [randn(2,N);ones(1,N)];
% Create a radial-basis function kernel
rbfKernel = @(x,y) exp(-norm(x-y)^2);
% Evaluate the kernel at all pairs of points
K = zeros(N,N);
HODLR_Mtrx = K;
for ii = 1:N
	for jj = 1:N
		K(ii,jj) = rbfKernel(dim2Points(1:2,ii),dim2Points(1:2,jj));
	end
end

% Plot the kernel matrix
%figure(1)
%imagesc(K)

% pick the first element in a random set to be the focus point.
focusPt = dim2Points(1:2,1);
% this point's distance will no longer be considered.
validPoints = dim2Points(1:2,:);
for ii = 1:N-1
  %scatter(dim2Points(1,:),dim2Points(2,:),[],dim2Points(3,:))
  % compute distances from the focus point to every other point in the set.

  % There is an issue on this line. distIdx allocates indices [1,99]; instead,
  % we need to permute the remaining indices of dim2Points.
  [~,distIdx] = sort(vecnorm(validPoints - focusPt));
  % sort the points in dim2Points into another set of 'valid points' we can
  % evaluate.
  validPoints = validPoints(:,distIdx);
  % evaluate the kernel against the focus point.
  rbfVector = zeros(1,length(validPoints));
  for jj = 1:length(validPoints)
    rbfVector(jj) = rbfKernel(validPoints(:,jj),focusPt);
  end
  HODLR_Mtrx(ii:N,ii) = rbfVector';
  validPoints = validPoints(:,2:end);
  focusPt = validPoints(:,1); % change the 1 to 'end'.
end

figure(2)
imagesc(HODLR_Mtrx)