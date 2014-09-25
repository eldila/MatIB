function Ux = Centered_Dx(U,dx)
% Centered_Dx  - Construct a centered difference approximation
% to the first derivative in the x direction with periodic BC.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

Ux = 0.5*(U(:,[2:end,1])-U(:,[end,1:end-1]))/dx;