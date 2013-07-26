% The purpose of this simulation is to test the decay rate of the maximum
% height of a fibre. We begin with a fibre wrapped around the periodic 
% domain with an initially sinusoidal profile.
% The fibre should oscillate with a decaying amplitude.  The rate of decay
% can be compute as in Chapter 3 of John Stockie's PhD thesis.  Here we
% assume that the resting length of the fibre is zero (L = 0).  
%

% Add PATH reference in order to run solver
addpath('../../solver/Peskin-TwoStep');
addpath('../../solver/utils');

% The number of grid points.
N = 64; 
Nb = 3*N;

% Parameter values.
L = 0;         % Resting length.
mu = 0.1;      % Viscosity.
sigma = 1;     % Spring constant.
rho = 1;       % Density.
A = 0.05;      % Initial height of the fibre.

% Time step and final time.
Tfinal = 1;
dt = 5e-3;
NTime = floor(Tfinal/dt);
Tfinal = dt*NTime;

% The initial velocity is zero.
IC_U = @(X,Y)zeros(size(X));
IC_V = @(X,Y)zeros(size(Y));

% The membrane is a sinusoidally perturbed line wrapping around the domain.
IC_ChiX = @(S)S;
IC_ChiY = @(S)0.5 + A*sin(2*pi*S);

% Get the maximum height of the fibre at each time.
Action = @( dx, dy, dt, indexT, X, Y, U, V, chiX, chiY, Fx, Fy)...
    PlotMembrane(X, Y, U, V, chiX,chiY,1);

% Do the IB solve.
[X, Y, S, U, V, chiX, chiY] = ...
    IBSolver(mu, rho, sigma, L, IC_U, IC_V, IC_ChiX, IC_ChiY, ...
    N, N, 1, 1, Nb, NTime, Tfinal, Action);

% Remove PATH reference to avoid clutter
rmpath('../../solver/Peskin-TwoStep');
rmpath('../../solver/utils');
