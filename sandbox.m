profile on
% make a kernel
rbf = @(e,r) exp(-(e*r).^2); ep = 5; d = 1;
% Pointervals at which to sample
x = CreatePoints(n^d,d,'u');
% compute the absolute difference
DM = DistanceMatrix(x,x);
% sample matrix
K_mtrx = rbf(ep,DM);
% option 1
n = 1000; l = 3; k = 10;
mask = MakeMask(n,l,k);
Omega = randn(n, 2*k);
Y = K_mtrx * (Omega .* mask);
profile off
profile viewer

% option 2
omega = randn(n, 2*k);
Y = Interlace(omega,n,l,k);

function mask = MakeMask(n,l,k)
    % build the first "row" of zeros and random values (Blue/white block in
    % Fig (12a)
    row = repmat([zeros(1,k) ones(1,k)], floor(n/(2^l)), 1);
    % Flip the existing block and stack the two to form a "block"
    row_flip = fliplr(row);
    block = [row; row_flip];
    % keep stacking blocks 2^(l-1) times to form the entire mask
    mask = repmat(block, 2^(l-1), 1);
end

% time turns flames to embers. You'll have new Septembers. In dreams, I meet
% you in warm conversation.
function output_mtrx = Interlace(omega,n,l,k)
end
