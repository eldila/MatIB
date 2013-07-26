function PlotMembrane(X, Y, U, V,chiX,chiY,figNum)
% Plots the membrane
%
% Assumption: Assuming chiX and chiY are column vectors
% Assumption: Assuming chiX(i+1)-chiX(i) < .5 and chiY(i+1)-chiY(i) < .5, for all points that don't cross the boundary
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

if nargin==7
    figure(figNum)
else
    figure;
end;

Lx = X(1,end)+X(1,2);
Ly = Y(end,1)+Y(2,1);

chiX = ShiftMembraneValues(chiX,Lx);
chiY = ShiftMembraneValues(chiY,Ly);

diffX = abs([chiX(2:end);chiX(1)]-chiX);
diffY = abs([chiY(2:end);chiY(1)]-chiY);

locX = find(diffX > Lx/2); % assuming x value can't change more than Lx/2 (without crossing boundary)
locY = find(diffY > Ly/2); % assuming y value can't change more than Ly/2 (without crossing boundary)
loc = sort(unique([locX;locY]));

clf;
axis([0 Lx 0 Ly]);
xlabel('x'); ylabel('y');

hold all;

quiver(X,Y,U,V);

if length(loc) > 0
   loc = [0;loc;length(chiX)];
   for i=2:length(loc)
       plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'r');
   end;
else
   plot([chiX;chiX(1)],[chiY;chiY(1)],'r');
end;

axis square;

drawnow;
hold off;

set(gca,'Box','on');
