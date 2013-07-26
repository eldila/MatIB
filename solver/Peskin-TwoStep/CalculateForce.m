function [Fx, Fy, fx, fy] = CalculateForce...
    (chiX, chiY, x, y, ds, dx, dy, Nb, Nx, Ny, Lx, Ly, sigma, L)
%
% This function computes the components of the forcing coming from
% interactions between the fluid and the immersed boundary.
% The components of the force are given by
% F(x,y) = int_Gamma f(s) * delta(x - chiX(s)) * delta(y - chiY(s)) * ds
% where 
% s parameterises the immersed boundary,
% f(s) is the force density on the boundary,
% chiX, chiY are the x and y positions of points on the membrane, 
% delta is a delta function.
%
% The force density is computed using a simple linear spring model with
% resting length L.  This leads to a force density of
%    f = [\chi_s * ( 1 - L / abs(\chi_s) ) ]_s
% where the subscript s refers to a derivative with respect to s, which
% parameterises the membrane.
%
% INPUTS:   chiX        A vector of length Nb that gives the current
%                       x-coordinates of the position of the immersed bounday
%                       at time-step n-1/2.
%           chiY        A vector of length Nb that gives the current
%                       y-coordinates of the position of the immersed bounday
%                       at time-step n-1/2.
%           x,y         Vectors of length N that give the x and y values
%                       occuring in the Eulerian grid.
%           ds          The step size in hte Lagrangian grid.
%           dx          The step size in x-direction for the Eulerian grid.
%           dy          The step size in y-direction for the Eulerian grid.
%           Nb          The total number of grid points in the Lagrangian
%                       grid.
%           Nx          The total number of grid points along x-dimension
%                       of the Eulerian grid.
%           Ny          The total number of grid points along y-dimension
%                       of the Eulerian grid.
%           Lx          The length of the domain along x-direction.
%           Ly          The length of the domain along y-direction.
%           sigma       The spring constant used to compute the force density.
%           L           The resting length of the membrane.
%
% OUTPUTS:  Fx, Fy      The x and y components of the forcing resulting from
%                       the fluid-membrane interactions.  These are N X N
%                       matrices that give the force at each point on the
%                       grid.
%            fx, fy     The x and y components of the force density on the
%                       membrane.  These are Nb X 1 vectors.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

% Compute the derivatives of the membrane positions with respect to s.
% These are centred differences evaluated at the points (s+1/2), which is
% equivalent to apply a forward difference to chi_s.
% For the forward difference on the periodic domain, the value of s_{Nb+1}
% is equivalent to the value of s_1.
dChiX = PeriodicDisplacement(chiX([2:Nb,1]), chiX, Lx) / ds;
dChiY = PeriodicDisplacement(chiY([2:Nb,1]), chiY, Ly) / ds;

% We also want the magnitude of the derivative.
dChi = sqrt(dChiX.^2 + dChiY.^2);

% Now we compute the tension, which is also given at the points (s+1/2),
% and points in the direction of the tangent vector.
Tx = sigma * dChiX .* (1 - L ./ dChi);
Ty = sigma * dChiY .* (1 - L ./ dChi);

% The force density is computed at the points s by taking a narrow centred
% difference of the tension. This is equivalent to computing a backward
% difference of T_{s+1/2}.
% This assumes a simple linear spring model with resting length L.
% For the backward difference on the periodic domain, the value of s_0
% is equivalent to the value of s_{Nb}.
fx = (Tx - Tx([Nb,1:Nb-1])) / ds;
fy = (Ty - Ty([Nb,1:Nb-1])) / ds;

% Convert the force density into a (sparse) diagonal matrices.
fxds = spdiags(fx * ds,0,Nb,Nb);
fyds = spdiags(fy * ds,0,Nb,Nb);

% Find out how many intervals of the grid the support of the 1D delta
% function will cover.
supp = delta1D();

% Get the indices of the points x_i, y_j where the 1D delta functions
% delta(x - chiX_l) and delta(y - chiY_l) are non-zero.
% These will be Nb X supp   and  supp X Nb  sized matrices, respectively.
indX = IndicesOfNonZeroDelta1D(chiX, Nb, Nx, dx, supp);
indY = IndicesOfNonZeroDelta1D(chiY, Nb, Ny, dy, supp)';

% We also need a matrix of indices that chi can take, which will need to
% have the same size as the other index matrices so we can do the
% subtraction at these points.
% This is Nb X supp.
indChi = repmat((1:Nb)',1,supp);

% Now we will construct sparse matrices to hold the values of the 1D delta
% functions at the points x_i - chiX_l
% where i = 1,..,Nx  and  l = 1,..,Nb
% stored in an Nb X Nx matrix 
% and at the points y_j - chiY_l
% where j = 1,..,Ny  and  l = 1,..,Nb
% stored in an Ny X Nb matrix 
deltaX = sparse(Nb, Nx);
deltaY = sparse(Ny, Nb);

% We only need to evaluate the distance between points in the grid and
% points on the membrane that will contribute to the delta function.  There
% is no point in computing all possible values of (x_i - chiX_l)  and 
% (y_j - chiY_l) when most of these will not contribute.
% The distances x_i - chiX_l are stored in an Nb X supp matrix.
% The distances y_j - chiY_l are stored in a supp X Nb matrix.
distX = PeriodicDistance(x(indX),chiX(indChi),Lx);
distY = PeriodicDistance(y(indY),chiY(indChi'),Ly);

% Now update the delta functions to include any non-zero contributions.
% These still have sizes (Nb X Nx) and (Ny X Nb).
% Each delta function is evaluated only at Nb * 4 points.  The other points
% will not contribute.
% These are still stored using a sparse data structure.
deltaX(sub2ind([Nb, Nx], indChi, indX)) = delta1D(distX, dx);
deltaY(sub2ind([Ny, Nb], indY, indChi')) = delta1D(distY, dy);

% Compute the components of the force via the (discrete) line integral
%  F(x,y) = int_Gamma f(s) * delta(x - chiX(s)) * delta(y - chiY(s)) * ds.
% We recall that each delta function (deltaX, deltaY) will have only Nb*4
% non-zero entires, while the force density matrices (fxds, fyds) are
% diagonal and have only Nb non-zero entries.
% The multiplication is performed using the sparse matrix structure, so
% this should be much less expensive than multiplying full matrices.
% The result will be an Ny X Nx matrix.
Fx = full(deltaY * fxds * deltaX);
Fy = full(deltaY * fyds * deltaX);

