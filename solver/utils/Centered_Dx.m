function D = Centered_Dx(Nx,Ny,Lx,Ly)
% Centered_Dx  - Construct a centered difference approximation
% to the first derivative in the x direction with periodic BC.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

dx = Lx/(Nx);

e = ones(Nx, 1);
Dx = spdiags([-1*e e],[-1 1], Nx, Nx);
Dx(1,end) = -1; Dx(end,1) = 1;
Dx = Dx/(2*dx);

I = eye(Ny);

D = kron(Dx,I); 
