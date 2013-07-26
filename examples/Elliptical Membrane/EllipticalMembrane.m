% In this example, we look at an elliptical membrane in a fluid that is
% initially at rest.  The ellipse should oscillate and eventually settle to
% a circle with a prescribed radius.

% Add PATH reference in order to run solver
addpath('../../solver/Peskin-TwoStep');
addpath('../../solver/utils');

% The number of grid points.
N = 2*round(2.^4); 
Nb = 3*N;

% Parameter values.
mu = 1;        % Viscosity.
sigma = 1e4;     % Spring constant.
rho = 1;       % Density.
rmin = 0.2;    % Length of semi-minor axis.
rmax = 0.4;    % Length of semi-major axis.
req = sqrt(rmin*rmax);
L = 0;  % Resting length.

% Time step and final time.
Tfinal = .04;
dt = 5e-5;
NTime = floor(Tfinal./dt)+1;
dt = Tfinal ./ NTime;

% The initial velocity is zero.
IC_U = @(X,Y)zeros(size(X));
IC_V = @(X,Y)zeros(size(Y));

% The membrane is an ellipse.
IC_ChiX = @(S)0.5 + rmax * cos(2*pi*S);
IC_ChiY = @(S)0.5 + rmin * sin(2*pi*S);

% The action function.
Action = @( dx, dy, dt, indexT, X, Y, U, V, chiX, chiY, Fx, Fy)...
    PlotMembrane(X, Y, U, V, chiX,chiY,1);

% Do the IB solve.
[X, Y, S, U, V, chiX, chiY] = ...
        IBSolver(mu, rho, sigma, L, IC_U, IC_V, IC_ChiX, IC_ChiY,...
        N, N, 1, 1, Nb, NTime, Tfinal, Action);
 
% Remove PATH reference to avoid clutter
rmpath('../../solver/Peskin-TwoStep');
rmpath('../../solver/utils');
