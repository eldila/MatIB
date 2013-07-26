function dist = PeriodicDistance(x, y, L)
% This function computes the distance between scalar values x and y as
% measured on the periodic interval [0,L].
% This means that the value x is equivalent to the point (x + n) for any
% integer n.
%
% INPUTS:  x, y      Matrices of the same size that give the values of x
%                    and corresponding values of y for which the distance
%                    is to be measured.  These values should all be between
%                    0 and L.
%             L      Length of the domain.
% 
% OUTPUTS: dist      A matrix (the same size as the inputs) that gives the
%                    distance between each x and y as measured on the
%                    periodic interval.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

% First compute difference between these.
dist = abs(x - y);

% Make sure we are measuring in the right direction.
dist = min(dist, L - dist);
