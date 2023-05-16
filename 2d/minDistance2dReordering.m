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
% Create a 2d C6 Matern kernel

%kernel = @(e,r) (1+e*r+2/5*(e*r).^2+1/15*(e*r).^3).*exp(-e*r); ep 10
%dim2Points = randn(n,2);
function [HODLR_Mtrx,focusPtIdx] = minDistance2dReordering(dim2Points,kMtrxFcn,ep)
HODLR_Mtrx = zeros(length(dim2Points));
focusPtIdx = zeros(length(dim2Points),1);
N = length(dim2Points);
% pick the first element in a random set to be the focus point.
focusPt = dim2Points(1,:);
focusPtIdx(1) = 1;
% this point's distance will no longer be considered.
validPoints = dim2Points(2:end,:);
    for ii = 2:N
      [~,distIdx] = sort(vecnorm(validPoints - focusPt,1,2));
      % sort the points in dim2Points into another set of 'valid points' we can
      % evaluate.
      validPoints = validPoints(distIdx,:);
      % evaluate the kernel against the focus point.
      DM = DistanceMatrix(validPoints,focusPt);
      rbfVector = kMtrxFcn(ep,DM);
      HODLR_Mtrx(ii:N,ii) = rbfVector;
      focusPt = validPoints(end,:); % change the 1 to 'end' for maximin.
      validPoints = validPoints(1:N-ii,:);
      focusPt(ii) = distIdx(end);
    end
end