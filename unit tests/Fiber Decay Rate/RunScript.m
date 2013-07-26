% The purpose of this simulation is to test the decay rate of the maximum
% height of a fibre. We begin with a fibre wrapped around the periodic 
% domain with an initially sinusoidal profile.
% The fibre should oscillate with a decaying amplitude.  The rate of decay
% can be compute as in Chapter 3 of John Stockie's PhD thesis.  Here we
% assume that the resting length of the fibre is zero (L = 0).  
%
% We provide some sample parameter values.
%
% Example 1:  mu = 0.1
%             sigma = 1
%             rho = 1
%             Decay rate (analytical) is -3.335
%             Frequency (analytical) is 8.583
%             Tfinal = 2
%
% Example 2:  mu = 0.1
%             sigma = 5
%             rho = 1
%             Decay rate (analytical) is -4.485
%             Frequency (analytical) is 21.305
%             Tfinal = 0.5
%
% Example 3:  mu = 0.5
%             sigma = 2
%             rho = 10
%             Decay rate (analytical) is -1.606
%             Frequency (analytical) is 3.754
%             Tfinal = 2
%

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
N = 2*round(2^5); 
Nb = 3*N;

% Parameter values.
mu = [0.1, 0.1, 0.5];      % Viscosity.
sigma = [1, 5, 2];     % Spring constant.
rho = [1, 1, 10];       % Density.
A = 0.05;      % Initial height of the fibre.

% Analytical decay rate and frequency. 
rate = [-3.335, -4.485, -1.606];
freq = [8.583, 21.305, 3.754];

% Time step and final time.
Tfinal = [2, .5, 2];
dt = 1e-4;
NTime = floor(Tfinal./dt)+1;
dt = Tfinal ./ NTime;

% Plot figure
h = figure;
t = (0:NTime(1)) * dt(1);
createfigure(h, t,A * exp(-1.1*t) .* abs(cos(6.1*t)));
saveas(h,'example.eps', 'epsc');

for i = 1:length(mu)
    % Run Simulation for Fibre Decay (for different spatial step-size)
    height = FibreDecayRates(A, mu(i), sigma(i), rho(i), NTime(i), Tfinal(i), N, Nb);

    % The maximum height obtained from the k = 2\pi mode (analytical).
    t = (0:NTime(i)) * dt(i);
    heightAnal = A * exp(rate(i) *t) .* abs(cos(freq(i)*t));

    % Plot the maximum height.
    h = figure;
    plot(t,[A,height],'b-');
    hold on;
    plot(t,heightAnal,'r--');
    hold off;
    xlabel('Time (s)');
    ylabel('Max Fibre Height (cm)');
    axis([0 NTime(i)*dt(i) 0 A]);
    legend('Computed','k = 2\pi Mode (Analytical)');
    saveas(h,sprintf('fig%d.eps', i), 'epsc'); % use epstopdf to convert to pdf    
    
    % Find the times where the fibre height is maximum.
    % These are the points where the derivative changes sign.
    tMax = find(and(height(2:end-1)-height(1:end-2)>0,...
        height(3:end)-height(2:end-1)<0)) + 1;

    % Get the maximum heights at those times.
    hMax = height(tMax);
    tMax = tMax * dt(i);

    % Estimate the decay rate and frequency if we have enough information.
    lambdaRe = [];
    lambdaIm = [];
    if length(tMax)>=2,
        lambdaRe = 1./(tMax(2:end)-tMax(1:end-1)) .* log(hMax(2:end)/hMax(1:end-1));
        lambdaIm = pi ./ (tMax(2:end)-tMax(1:end-1));
    end;

    fprintf('Case %d: mu=%f, sigma=%f, rho=%f \n', i, mu(i), sigma(i), rho(i)) ; 
    fprintf('Decay Rate\n'); 
    fprintf('Computation: %f, Asymptotic:%f \n', lambdaRe(1), rate(i)) ; 
    fprintf('Decay Frequency\n'); 
    fprintf('Computation: %f, Asymptotic:%f \n\n', lambdaIm(1), freq(i)) ; 
    lambdaRe
    lambdaIm
    
end;



% Remove PATH reference to avoid clutter
rmpath('../../solver/Peskin-TwoStep');
rmpath('../../solver/utils');