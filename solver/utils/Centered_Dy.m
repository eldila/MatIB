function D = Centered_Dy(Nx, Ny)
% Centered_Dy  - Construct a centered difference approximation
% to the first derivative in the y direction with periodic BC.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

dy = 1/(Ny);

e = ones(Ny, 1);
Dy = spdiags([-1*e e],[-1 1], Ny, Ny);
Dy(1,end) = -1; Dy(end,1) = 1;
Dy = Dy/(2*dy);

I = eye(Nx);

D = kron(I,Dy); 