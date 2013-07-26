function disp = PeriodicDisplacement(x, y, L)
% This function computes the displacement between to values x and y which
% lie on the periodic interval [0,L].
% Here we assume that the absolute value of the displacement can never be
% larger than L/2 - if it is, we are measuring in the wrong direction.
% For example, when L=1, the value of 
%   x - y = 0 - 0.9
% should be equal to 0.1
% since we suppose that that value of zero (which is equivalent to 1) is to
% the right of 0.9.
% Similarily, the value of 
%   x - y = 0.9 - 0
% should be equal to -0.1.
%
% INPUTS:   x, y    Matrices of the same size that give the points for
%                   which the discplacement is required.  These values
%                   should all be between 0 and 1.
%              L    The length of the interval.
%
% OUTPUTS:  disp    A matrix the same size as x and y, which gives the
%                   displacement on the periodic interval.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

% Compute the displacement.
% If the value of x - y is greater than L/2, we suppose that y needs to be
% shifted so that it is to the right of x; this amounts to increasing the
% value of y by L.
% If the value of x - y is less than -L/2, we suppose that y needs to be
% shifted so that it is to the left of x; this amounts to decreasing the
% value of y by L.
disp = (x - y) + L*(x - y < -L/2) - L*(x - y > L/2);
