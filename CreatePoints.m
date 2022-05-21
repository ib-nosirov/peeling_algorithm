function [points, N] = CreatePoints(N,d,gridtype)
% Computes a set of N points in [0,1]^d
% Note: could add variable interval later
% Inputs:
% N: number of interpolation points
% d: space dimension
% gridtype: 'c'=Chebyshev, 'f'=fence(rank-1 lattice),
%    'h'=Halton, 'l'=latin hypercube, 'r'=random uniform, 
%    's'=Sobol, 'u'=uniform grid
% Outputs:
% points: an Nxd matrix (each row contains one d-D point)
% N: might be slightly less than original N for 
%    Chebyshev and gridded uniform points
% Calls on: chebsamp, lattice, lhsamp, gridsamp
% Also needs: fdnodes, gaussj
% Requires Statistics Toolbox for haltonset and sobolset.

switch gridtype
    case 'c'
        ppd = zeros(1,d);
        for j=1:d
            ppd(j) = floor(nthroot(N,d+1-j));
            N = N/ppd(j);
        end
        gam = 0.5*ones(1,d);  % density for point distribution, 0.5=Chebyshev
        points = chebsamp([zeros(1,d); ones(1,d)], ppd, gam);
        N = prod(ppd);
    case 'f'
        points = lattice(N,d);  % N should be(?) power of 2
    case 'h' 
        temp = haltonset(d);
        points = net(temp,N);
    case 'l'
        points = lhsamp(N,d);
    case 'r'
        rand('state',47); 
        points = rand(N,d);
    case 's'
        temp = sobolset(d);
        points = net(temp,N);
        points = N*points/(N-1);        
    case 'u'
        ppd = zeros(1,d);
        for j=1:d
            ppd(j) = floor(nthroot(N,d+1-j));
            N = N/ppd(j);
        end
        points = gridsamp([zeros(1,d); ones(1,d)], ppd);
        N = prod(ppd);
    otherwise
        error('Please use c, f, h, r, s or u data types')
end
