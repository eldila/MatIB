function Uyy = Centered_Dyy(U,dy)
% Centered_Dyy  - Construct a centered difference approximation
% to the second derivative in the y direction with periodic BC.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

Uyy = (U([2:end,1],:)+U([end,1:end-1],:)-2*U)/dy^2;