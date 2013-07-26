function delta = delta1D(x,dx)
% 
% This function computes a discrete approximation to the 1D Dirac-delta
% function.
% The discrete approximation should have a support contained in
% [x-dx,x+dx].
% The output will be have the same size as the input x.
% If no inputs are provided, this function will instead return an integer
% describing the number of grid points covered by the support of this
% function.  This integer (r) will be even, and the support will then be
% contained in the interval [-r/2,r/2].
%
% INPUTS:    x      The values at which the delta function is to be
%                   evaluated.  This can be a scalar, vector, or matrix.
%            dx     The step size of the grid, which also provides
%                   information about the support of the discrete delta
%                   function.
%
% OUTPUTS:   delta  A matrix the same size as the input x, which gives
%                   the values of the discrete delta function on the grid.
%                   If no inputs are provided, this will instead return an
%                   integer describing the size of the support of the
%                   discrete delta function.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

if nargin < 2
    
    % Describe the support.
    delta = 4;
    
else
    
    % Compute the delta function.
    r = abs(x)/dx;
    delta = ( (3 - 2*r + sqrt(1 + 4*r - 4*r.*r) ) / (8*dx) ) .* (r < 1) ;
    delta = ( (5 - 2*r - sqrt(-7 + 12*r - 4*r.*r) ) / (8*dx) ) .* (r < 2).*(r >= 1) + delta;
    
end;
