% Copyright (c) 2020, Yaqing Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

% Define surface type: 1 - hexagonal
surface_type = 1;

% Load the user parameters and jump vectors
if surface_type == 1
    hex_ui
    hex_vectors
else
    % for future support
    % in additional surfaces
end

surf_gen(lattice,Temp,tau,E);

plot_isf(Azmth,azmthStr,dK)

