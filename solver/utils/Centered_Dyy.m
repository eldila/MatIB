function D = Centered_Dyy(Nx, Ny)
% Centered_Dyy  - Construct a centered difference approximation
% to the second derivative in the y direction with periodic BC.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

dy = 1/(Ny);

e = ones(Ny, 1);
Dyy = spdiags([e -2*e e],[-1 0 1], Ny, Ny);
Dyy(1,end) = 1; Dyy(end,1) = 1;
Dyy = Dyy/(dy*dy);

I = eye(Nx);

D = kron(I,Dyy); 