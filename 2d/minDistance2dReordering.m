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

function HODLR_Mtrx = minDistance2dReordering(dim2Points,kMtrxFcn)
HODLR_Mtrx = zeros(length(dim2Points));
N = length(dim2Points);
% pick the first element in a random set to be the focus point.
focusPt = dim2Points(1,:);
% this point's distance will no longer be considered.
validPoints = dim2Points(2:end,:);
    for ii = 2:N
      [~,distIdx] = sort(vecnorm(validPoints - focusPt,1,2));
      % sort the points in dim2Points into another set of 'valid points' we can
      % evaluate.
      validPoints = validPoints(distIdx,:);
      % evaluate the kernel against the focus point.
      rbfVector = zeros(size(validPoints,1),1);
      for jj = 1:size(validPoints,1)
          rbfVector(jj) = kMtrxFcn(validPoints(jj,:),focusPt);
      end
      HODLR_Mtrx(ii:N,ii) = rbfVector;
      focusPt = validPoints(end,:); % change the 1 to 'end' for maximin.
      validPoints = validPoints(1:N-ii,:);
    end
end