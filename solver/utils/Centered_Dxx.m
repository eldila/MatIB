function D = Centered_Dxx(Nx, Ny, Lx, Ly)
% Centered_Dxx  - Construct a centered difference approximation
% to the second derivative in the x direction with periodic BC.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

dx = Lx/(Nx);

e = ones(Nx, 1);
Dxx = spdiags([e -2*e e],[-1 0 1], Nx, Nx);
Dxx(1,end) = 1; Dxx(end,1) = 1;
Dxx = Dxx/(dx*dx);

I = eye(Ny);

D = kron(Dxx,I); 