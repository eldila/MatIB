function chi = ShiftMembraneValues(chi, L)
%
% This shifts the values of a membrane position so that they lie in the
% interval [0,L].  Here we are assuming that the domain is periodic.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

% Do the shift.
chi = mod(chi, L);
