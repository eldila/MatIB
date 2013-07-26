function rAv = EllipticalMembrane(rmin, rmax, mu, sigma, rho, L, NTime, Tfinal, N, Nb)
% Solves the elliptical membrane problem for the inputted parameters
% See documentation or RunScript.m to demonstrate the use of the function.


% The initial velocity is zero.
IC_U = @(X,Y)zeros(size(X));
IC_V = @(X,Y)zeros(size(Y));

% The membrane is an ellipse.
IC_ChiX = @(S)0.5 + rmax * cos(2*pi*S);
IC_ChiY = @(S)0.5 + rmin * sin(2*pi*S);

% The action function.
% Get the average radius of the ellipse at each time.
global radius
radius = zeros(1,NTime);
Action = @( dx, dy, dt, indexT, X, Y, U, V, chiX, chiY, Fx, Fy)...
    GetRadius(indexT, chiX, chiY);
rAv = zeros(length(N),NTime);

% See how the radius varies with time.
for jj = 1:length(N)
    
    % Do the IB solve.
    [X, Y, S, U, V, chiX, chiY] = ...
        IBSolver(mu, rho, sigma, L, IC_U, IC_V, IC_ChiX, IC_ChiY,...
        N(jj),N(jj), 1, 1, Nb(jj), NTime, Tfinal, Action);
    
    % Store the radius estimates from this value of N.
    rAv(jj,:) = radius;
    radius(end)
end;

function GetRadius(indexT, chiX, chiY)
% This function determines the radius of a circular membrane at the time 
% step indicated by indexT.
% The result is stored in the corresponding location in the global variable
% radius.

% The radius variable is global.
global radius

% Compute the maximum height at the current time step.
radius(indexT) = mean(sqrt((chiX-mean(chiX)).^2 + (chiY-mean(chiY)).^2));
