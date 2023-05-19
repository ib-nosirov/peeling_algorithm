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

function [HODLR_Mtrx,focusPtIdx] = minDistance2dReordering(dim2Points,kMtrxFcn,ep,HODLR_Mtrx)
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
      focusPt = validPoints(end,:); % change the 1 to 'end' for maximin.
      validPoints = validPoints(1:N-ii,:);
      focusPtIdx(ii) = distIdx(end);
    end
    HODLR_Mtrx = HODLR_Mtrx(focusPtIdx,focusPtIdx); % check that focusPtIdx has all of the indices.
    % pick points that are far away from all points picked so far; current
    % approach is too greedy.

    % things to try: create a constellation of 5-6 points and work it out
    % by hand.
end