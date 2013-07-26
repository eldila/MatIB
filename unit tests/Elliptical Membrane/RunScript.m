% In this example, we look at an elliptical membrane in a fluid that is
% initially at rest.  The ellipse should oscillate and eventually settle to
% a circle with a prescribed radius.

% Add PATH reference in order to run solver
addpath('../../solver/Peskin-TwoStep');
addpath('../../solver/utils');

% Set Figure settings
set(0,'defaultaxesfontsize',20);
set(0,'defaultaxeslinewidth',0.75);
set(0,'defaultlinelinewidth',2);
set(0,'defaultpatchlinewidth',1);
set(0,'defaultlinemarkersize',10);

% The number of grid points.
N = 2*round(2.^(5:7)); 
Nb = 4*N;

% Parameter values.
mu = 0.05;     % Viscosity.
sigma = 0.41;  % Spring constant.
rho = 1;       % Density.
rmin = .35;    % Length of semi-minor axis.
rmax = .25;    % Length of semi-major axis.
L = 0.0;       % Resting length.

% Time step and final time.
Tfinal = 1.0;
dt = 0.0001;
NTime = floor(Tfinal./dt)+1;

% Run Simulation for Elliptical Membrane
req = sqrt(rmin*rmax);
rAv = EllipticalMembrane(rmin, rmax, mu, sigma, rho, L, NTime, Tfinal, N, Nb);

% Plot the average radius.
colour = [' b-', 'r--', 'g-.', ' k:'];
t = (1:NTime)*dt;

leg = [];
h = figure;
fig = zeros(length(N)+1,1);
for jj = 1:length(N)
   fig(jj) = plot(t,rAv(jj,:),colour(mod((jj-1)*3+(0:2),12)+1));
    hold on;
    leg{end+1} = ['N = ', num2str(N(jj))];
end;
plot(t,req*ones(size(t)),colour(mod((length(N))*3+(0:2),12)+1));
leg{end+1} = 'Equilibrium';
hold off;
xlabel('Time (s)');
ylabel('Average Radius (cm)');
leg = legend(leg);
saveas(h,'fig.eps', 'epsc'); % use epstopdf to convert to pdf  

% Remove PATH reference to avoid clutter
rmpath('../../solver/Peskin-TwoStep');
rmpath('../../solver/utils');