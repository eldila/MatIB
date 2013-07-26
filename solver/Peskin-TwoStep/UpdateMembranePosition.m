 function [chiX, chiY] = UpdateMembranePosition...
     (U, V, chiX, chiY, chiXConv, chiYConv, x, y, dt, dx, dy, Nb, Nx, Ny, Lx, Ly)
% 
% This function updates the position of the immersed boundary by evolving
% the equation 
%     Chi_t = int_\Omega u(x,t) delta(x - Chi(s,t)) dx.
%
% INPUTS:  U,V                 The x and y components of the fluid velocity at the
%                              previous time step, in matrix format.
%          chiX, chiY          The x and y component of the membrane position at the
%                              previous time step n-1.  These are vectors.
%          chiXConv, chiYConv  The x and y component of the membrane position
%                              used in the integral (summation). These are vectors.                     
%          x, y                The x and y positions in the grid.  These are vectors
%                              of length N.
%          dt                  The time step, a scalar.
%          dx                  The grid spacing, a scalar, in x-direction.
%          dy                  The grid spacing, a scalar, in y-direction.
%          Nb                  The total number of points used to describe the
%                              membrane position.
%          Nx                  The total number of grid points along x-direction
%                              of the domain.
%          Ny                  The total number of grid points along y-direction
%                              of the domain.
%          Lx                  The length of the domain in x-direction.
%          Ly                  The length of the domain in y-direction.
%
% OUTPUTS: chiX, chiY The updated x and y components of the membrane
%                     position.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

% The length of the support of the discrete delta function.  
supp = delta1D();

% Get the indices of the points x = ii*dx, y = jj*dx where the delta
% function is non-zero.  These are Nb X (supp^2) matrices.
[ii,jj] = IndicesOfNonzeroDelta(chiXConv, chiYConv, Nb, Nx, Ny, dx, dy, supp);

% Get the corresponding index to locate these points.
ind = sub2ind([Ny,Nx],jj,ii);

% We need to shift the values of chi into the square and augment the
% vectors to make them the same size as x and y for the subtraction.
% These will be matrices of size Nb X (supp^2).
chiXConvBig = repmat(ShiftMembraneValues(chiXConv,Lx),1,supp^2);
chiYConvBig = repmat(ShiftMembraneValues(chiYConv,Ly),1,supp^2);

% Compute the 1D delta functions.
deltaX = delta1D(PeriodicDistance(x(ii), chiXConvBig, Lx), dx);
deltaY = delta1D(PeriodicDistance(y(jj), chiYConvBig, Ly), dx);

% Update the membrane position.
chiX = chiX + (dt*dx*dy) * sum(U(ind) .* deltaX .* deltaY, 2);
chiY = chiY + (dt*dx*dy) * sum(V(ind) .* deltaX .* deltaY, 2);

% Now correct so that all values are in [0,1).
chiX = ShiftMembraneValues(chiX, Lx);
chiY = ShiftMembraneValues(chiY, Ly);

