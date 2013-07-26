function height = FibreDecayRates(A, mu, sigma, rho, NTime, Tfinal, N, Nb)
% Solves the Fibre Decay problem for the inputted parameters
% See documentation or RunScript.m to demonstrate the use of the function.

% Parameters
L = 0;         % Resting length.

% The initial velocity is zero.
IC_U = @(X,Y)zeros(size(X));
IC_V = @(X,Y)zeros(size(Y));

% The membrane is a sinusoidally perturbed line wrapping around the domain.
IC_ChiX = @(S)S;
IC_ChiY = @(S)0.5 + A*sin(2*pi*S);

% We want to store the maximum fibre height.
global glob_height;
glob_height = zeros(1,NTime);

% Get the maximum height of the fibre at each time.
Action = @( dx, dy, dt, indexT, X, Y, U, V, chiX, chiY, Fx, Fy)...
   GetFibreHeight(indexT, chiY);

% Do the IB solve.
[X, Y, S, U, V, chiX, chiY] = ...
    IBSolver(mu, rho, sigma, L, IC_U, IC_V, IC_ChiX, IC_ChiY,...
             N, N, 1, 1, Nb, NTime, Tfinal, Action);
    
height = glob_height;


function GetFibreHeight(indexT, chiY)
% This function determines the maximum height of a fibre (as measured from
% the position chiY = 0.5) at the time step indicated by indexT.
% The result is stored in the corresponding location in the global variable
% height.
% This function accompanies the simulation ComputeDecayRates.m

% The height variable is global.
global glob_height;

% Compute the maximum height at the current time step.
glob_height(indexT) = max(abs(chiY-0.5));
