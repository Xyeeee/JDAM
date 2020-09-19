% Copyright (c) 2020,Yaqing Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

function surf_gen(lattice,temperature,tau,E)

%%  Params.
SE_k1=9.6485335; % meV/Angstrm to a.m.u * Angstrm/(ps)^2 ... 1.6021766208e-19 / 1000 * 1/1.660539e-27 * 1e20 / 1e24
SE_kB_nano = 0.8317035; % Boltzmann constant in Aˆ2 amu psˆ-2 Kˆ-1 .... 8.6173303e-5 eV/K * k1
Boltzmann_const = SE_kB_nano/SE_k1; % In meV/K
system.T = temperature; % Unit: Kelvin

%% Reduce dimentions for non-activated sites
% If Tau contains Infs or NaNs, reduce the dimention of the problem
% accordingly
dummyTau = triu(tau)+tril(tau');
tmp = isinf(dummyTau) | isnan(dummyTau);
Col = find(~all(tmp));
Row = find(~all(tmp'));

tau = tau(Row,Col);
lattice.m = size(tau,1);
lattice.singles = lattice.singles(Row,Col);
lattice.s = lattice.s(Row,Col,:,:);
E = E(Row);

%% Calculate Lamda, concentraions and completion of Tau

% Calculate energy differences between sites
lattice.E_diff = zeros(lattice.m); %DEFAULT Site Energies
for i = 1:lattice.m
    for j = 1:lattice.m
        
        % the energy difference between site i and site j
        lattice.E_diff(i,j) = E(j)-E(i);
        
    end
end

% Calculate the lambda matrix
Gamma = @(e) exp(-e/(Boltzmann_const*system.T)); 
lattice.lambda = Gamma(lattice.E_diff); % The ratio of site(energy state) concentrations given thermodynamically

% Calculate the relative site concentrations
lattice.c = lattice.lambda(1,:)./sum(lattice.lambda(1,:));%Normalized site concentration vectors

% Calculate the lower triangular part of tau matrix
% tau as input is an mxm upper triangular matrix.
% Theory: tau(j,i) = tau(i,j)*lambda(i,j)
lattice.tau = tau;
for i = 1:lattice.m
    for j = (i+1):lattice.m
        lattice.tau(j,i) = lattice.tau(i,j)*lattice.lambda(i,j);
    end
end

%% Calculating each individual jump rate

jump_rates = zeros(lattice.m,lattice.m,max(max(lattice.singles)));
for i = 1:lattice.m
    for j = 1:lattice.m
        for k = 1:lattice.singles(i,j)
            jump_rates(i,j,k) = 1/(lattice.tau(i,j)*lattice.singles(i,j));
        end
    end
end

%% Saving the surface configurations for plotting
save('surface');
end
