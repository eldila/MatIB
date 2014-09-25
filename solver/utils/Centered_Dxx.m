function Uxx = Centered_Dxx(U,dx)
% Centered_Dxx  - Construct a centered difference approximation
% to the second derivative in the x direction with periodic BC.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

Uxx = (U(:,[2:end,1])+U(:,[end,1:end-1])-2*U)/dx^2;