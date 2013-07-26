function ii = IndicesOfNonZeroDelta1D(chi, Nb, N, dx, supp)
%
% This function computes the indices ii of the points x = ii*dx 
% where a 1D discrete delta function of the form
% delta(x-chiX) may be nonzero.
% Here the support of the 1D delta functions is contained in an interval of
% length (supp * dx), centred at the origin.
%
% INPUTS:  chi         Vectors of length Nb that indicate one coordinate of
%                      the membrane position.
%          Nb          The number of grid points on the membrane.
%          N           The number of grid points along each dimension of
%                      the domain.
%          dx          The spatial resolution of the Eularian grid.
%          supp        An integer describing the size of the support of the
%                      discrete 1D delta function.  This should be even.
%
% OUTPUTS: ii          A matrix of size Nb X supp that givex the x 
%                      indices where the 1D discrete delta function is
%                      non-zero.
%
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

% Make sure the support is an even number.
supp = ceil(supp/2)*2;

% Get the index of the grid point to the lower left of each membrane point.
% This will be an Nb X 1 vectors.
ii = floor(chi/dx + 1);

% Get all the different x indices that must be considered.
% This will be an Nb X supp matrix.
ii = repmat(ii,1,supp) + repmat(-supp/2+1:supp/2,Nb,1);

% Make sure all indices are between 1 and N.
ii = mod(ii-1,N) + 1;
