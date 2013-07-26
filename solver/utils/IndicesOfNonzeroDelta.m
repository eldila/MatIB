function [ii,jj] = IndicesOfNonzeroDelta(chiX, chiY, Nb, Nx, Ny, dx, dy, supp)
%
% This function computes the indices (ii,jj) of the points x = ii*dx, 
% y = jj*dx where a 2D discrete delta function of the form
% delta(x-chiX)*delta(y-chiY) may be nonzero.
% Here the support of the 1D delta functions is contained in an interval of
% length (supp * dx), centred at the origin.
%
% INPUTS:  chiX,chiY   Vectors of length Nb that indicate the x and y
%                      coordinates of the membrane position.
%          Nb          The number of grid points on the membrane.
%          Nx          The number of grid points along x-direction of
%                      the domain.
%          Ny          The number of grid points along y-direction of
%                      the domain.
%          dx          The spatial resolution of the Eularian grid 
%                      along x-direction.
%          dy          The spatial resolution of the Eularian grid 
%                      along y-direction.
%          supp        An integer describing the size of the support of the
%                      discrete 1D delta function.  This should be even.
%
% OUTPUTS: ii, jj      Matrices of size Nb X (supp^2) that give the x and y
%                      indices where the 2D discrete delta function is
%                      non-zero.
%
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

% Make sure the support is an even number.
supp = ceil(supp/2)*2;

% Get the index of the grid point to the lower left of each membrane point.
% These will be Nb X 1 vectors.
ii = floor(chiX/dx + 1);
jj = floor(chiY/dy + 1);

% Get all the different x indices that must be considered.
% This will be an Nb X (supp^2) matrix.
% Each row lists the required x-indices, which are then repeated (supp)
% times since we will need to consider all different combinations and the x
% and y values.
% For example, if supp = 4 and the index to the immediate left of the first
% membrane point is (i0 = 2), then the first row of ii will be
% (1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4).
ii = repmat(ii,1,supp) + repmat(-supp/2+1:supp/2,Nb,1);
ii = repmat(ii,1,supp);

% Get all the different y indices that must be considered.
% This will be an Nb X (supp^2) matrix.
% Each row lists the required y-indices, each of which is repeated (supp)
% times to enable us to consider all different combinations of x and y.
% For example, if supp = 4 and the index immediately below the first
% membrane point is (j0 = 2), then the first row of jj will be
% (1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4).
modJ = reshape(repmat((-supp/2+1:supp/2)',1,supp*Nb)',Nb,supp^2);
jj = repmat(jj,1,supp^2) + modJ;

% Make sure all indices are between 1 and N.
ii = mod(ii-1,Nx) + 1;
jj = mod(jj-1,Ny) + 1;
