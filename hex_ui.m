% Copyright (c) 2020, Yaqing Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

% ============ User interface - for a hexagonal surface ============
% Computation based on the paper by Tuddenham et. al.,Surface Science 604 (2010) 1459–1475
% Note: In the current version only single jumps are considered.


%% Surface azimuths for calc and plot

dK = linspace(0,2.5); % Parallel momentum transfer

Azmth = [1,0;0,1]; % surface azimuths to be considered, Azmth(i,:) is the i'th azimuth

% Legends for the azimuths
azmthStr = {'$[11\bar{2}]$','$[1\bar{1}0]$'};

%% Temperature and Site energies [K,meV]
% This will affect the site concentrations determined thermodynamically,
% i.e. the ratio of probabilities of residence for a particle on two
% different sites is proportional to exp(-Energy_differece/kT) where k is
% the Boltzmann constant and T is the system temperature in Kelvin

Temp = 300;      % [K]

fcc   = 46;      % [meV]
hcp   = 0 ;      % [meV]
bri_1 = 300; % [meV]
bri_2 = bri_1; % [meV]
bri_3 = bri_1; % [meV]
top   = 300; % [meV]

E = [top,bri_1,bri_2,bri_3,hcp,fcc];

%% tau matrix [picoseconds]
% sites in one unit cell has dimension 6 by 6. The total jump rate from a
% site of type i to a destination site of type j, both indices go from 1 to 6
% e.g. fcc to hcp, is given by tau(6,5)^(-1)

% NOTE: Only the UPPER TRIANGULAR part of tau needs to be initialized
% because tau(j,i) = tau(i,j)*c(i)/c(j) where c is the normalized site 
% particle concentration vector

% An example tau may be initialized with the same magnitude in unit picosecond
tau = ones(6)*3; %ps

% ===========================================
% To limit jumps to certain sites, you may configure the taus to be all
% infinities apart from the sites allowed. , for example in case of hollow
% sites only

% % Include only top sites
tau = ones(6).*inf; % setting all elements to be infinity
tau(1,1) = 3;

% % Include only hollow sites
% tau = ones(6).*inf; % setting all elements to be infinity
% tau(5,6) = 10;
% 
% 
% % Include only bridge sites
% tau = ones(6).*inf;
% tau(2,3) = 1;
% tau(2,4) = 1;
% tau(3,4) = 1;

% % Include only top and bridge sites
% tau = ones(6).*inf;
% tau(1:4,1:4) = 1;


% ===========================================


