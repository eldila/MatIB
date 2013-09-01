function Uy = Centered_Dy(U,dy)
% Centered_Dy  - Construct a centered difference approximation
% to the first derivative in the y direction with periodic BC.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

Uy = 0.5*(U([2:end,1],:)-U([end,1:end-1],:))/dy;