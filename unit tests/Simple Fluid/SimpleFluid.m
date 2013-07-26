function [SupErrorU, SupErrorV, L2ErrorU, L2ErrorV] = SimpleFluid(mu, rho, L, NTime, Tfinal, N)
% Solves a simple fluid problem for the inputted parameters
% See documentation or RunScript.m to demonstrate the use of the function.


% The initial velocity is zero.
IC_U = @(X,Y) 1 - 2*cos(2*pi*X).*sin(2*pi*Y);
IC_V = @(X,Y) 1 + 2*sin(2*pi*X).*cos(2*pi*Y);

% The membrane is an ellipse. 
IC_ChiX = @(S)0.5 + .25 * cos(2*pi*S);
IC_ChiY = @(S)0.5 + .25 * sin(2*pi*S);
Nb = 3*N;

% Error Vectors
SupErrorU = zeros(size(N));
SupErrorV = zeros(size(N));
L2ErrorU = zeros(size(N));
L2ErrorV = zeros(size(N));

% Action Function
Action = @( dx, dy, dt, indexT, X, Y, U, V, chiX, chiY, Fx, Fy) dx;

% See how the radius varies with time.
for jj = 1:length(N)
    
    % Do the IB solve.
    [X, Y, S, U, V, chiX, chiY] = ...
        IBSolver(mu, rho, 0.0, 0.0, IC_U, IC_V, IC_ChiX, IC_ChiY,...
        N(jj),N(jj), 1, 1, Nb(jj), NTime, Tfinal, Action);
    
    % Exact Solution
    UExact = 1 - 2*exp(-8*pi*pi*Tfinal*mu/rho)*cos(2*pi*(X-Tfinal)).*sin(2*pi*(Y-Tfinal));
    VExact = 1 + 2*exp(-8*pi*pi*Tfinal*mu/rho)*sin(2*pi*(X-Tfinal)).*cos(2*pi*(Y-Tfinal));

    % Norm Errors
    SupErrorU(jj) = max(abs(U(:) - UExact(:)));
    SupErrorV(jj) = max(abs(V(:) - VExact(:)));
    L2ErrorU(jj) = norm(U(:) - UExact(:),2)/N(jj);
    L2ErrorV(jj) = norm(V(:) - VExact(:),2)/N(jj);
end;

