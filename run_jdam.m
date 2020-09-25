% Copyright (c) 2020, Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

%% Clearing previous relavant workspace 

clear params;
clear lattice;
%% Define Parameters

params.tSE = 0:600;

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

[lattice,params] = surf_gen(lattice,params);

[Wp_all,D_abs,ISF] = plot_isf(lattice,params);


