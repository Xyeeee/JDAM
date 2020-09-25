% Copyright (c) 2020, Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

% ============ User interface - for a hexagonal surface ============
% Computation based on the paper by Tuddenham et. al.,Surface Science 604 (2010) 1459–1475
% Note: The Tuddenham paper use scaling factors for representability; those
% are implemented here as and option.
% Note: In the current version only single jumps are considered.

global scaleDecayInPlot

scaleDecayInPlot = 1;

%% Surface azimuths for calc and plot

lattice.a = 1;%2.71; % lattice constant
params.dKscale = 2*pi;%1; % a scaling factor for presentability; it multiply the dK parameter, then is being devided by at plotting

params.dK = [0:0.01:2.5]*params.dKscale; % Parallel momentum transfer

params.Azmth = [1,0 ; 0,1]; % surface azimuths to be considered, Azmth(i,:) is the i'th azimuth

% Legends for the azimuths
params.azmthStr = {'$[11\bar{2}]$','$[1\bar{1}0]$'}; % Two high symmetry directions for hexagonal surfaces

%% 
params.Psurvival = -1; % Survival Probability: if set to other than -1, 
% must also set other parameters to be explicitly consistent with jumps allowed on hollow sites only.

%% Temperature and Site energies [K,meV]
% This will affect the site concentrations determined thermodynamically,
% i.e. the ratio of probabilities of residence for a particle on two
% different sites is proportional to exp(-Energy_differece/(kB*T)) where kB is
% the Boltzmann constant and T is the system temperature in Kelvin
% Boltzmann constant in units of meV/K is defined in the file surf_gen.m

params.T = 300;      % [K]

fcc   = 46;      % [meV]
hcp   = 0 ;      % [meV]
bri_1 = 300; % [meV]
bri_2 = bri_1; % [meV]
bri_3 = bri_1; % [meV]
top   = 300; % [meV]

params.E = [top,bri_1,bri_2,bri_3,hcp,fcc];

%% tau matrix [picoseconds]
% sites in one unit cell has dimension 6 by 6. The total jump rate from a
% site of type i to a destination site of type j, both indices go from 1 to 6
% e.g. fcc to hcp, is given by tau(6,5)^(-1)

% NOTE: Only the UPPER TRIANGULAR part of tau needs to be initialized
% because tau(j,i) = tau(i,j)*c(i)/c(j) where c is the normalized site 
% particle concentration vector

% An example tau may be initialized with all ones in unit picosecond
params.tau = ones(6)*1; %ps

% ===========================================
% To limit jumps to certain sites, you may configure the taus to be all
% infinities apart from the sites allowed. For example in case of hollow
% sites only

% Include only top sites
params.tau = ones(6).*inf; % setting all elements to be infinity
params.tau(1,1) = 3;

% % Include only hollow sites
% params.tau = ones(6).*inf; % setting all elements to be infinity
% params.lambda = [1 1/9.5; 9.5 1]';
% params.tau(5,6) = 20;
% params.Psurvival = 0.7; % Survival Probability: if set to other than -1, must also set other parameters to confrom to the explicit jump on hollow sites case.

 
% % Include only bridge sites
% params.tau = ones(6).*inf;
% params.tau(2,3) = 1;
% params.tau(2,4) = 1;
% params.tau(3,4) = 1;

% % Include only top and bridge sites
% params.tau = ones(6).*inf;
% params.tau(1:4,1:4) = 1;


% ===========================================


