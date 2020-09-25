% Copyright (c) 2020, Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

function [lattice,params] = surf_gen(lattice,params)

%%  Params.
SE_k1=9.6485335; % meV/Angstrm to a.m.u * Angstrm/(ps)^2 ... 1.6021766208e-19 / 1000 * 1/1.660539e-27 * 1e20 / 1e24
SE_kB_nano = 0.8317035; % Boltzmann constant in Aˆ2 amu psˆ-2 Kˆ-1 .... 8.6173303e-5 eV/K * k1
params.kB = SE_kB_nano/SE_k1; % In meV/K

%% Reduce dimentions for non-activated sites
% If Tau contains Infs or NaNs, reduce the dimention of the problem
% accordingly
dummyTau = triu(params.tau)+tril(params.tau');
tmp = isinf(dummyTau) | isnan(dummyTau);
Col = find(~all(tmp));
Row = find(~all(tmp'));

params.tau = params.tau(Row,Col);
lattice.m = size(params.tau,1);
lattice.singles = lattice.singles(Row,Col);
lattice.s = lattice.s(Row,Col,:,:);
params.E = params.E(Row);

%% Calculate Lamda (if required), and calculate concentraions and completion of Tau
if isfield(params,'E') && ~isfield(params,'lambda')
    
    % Calculate energy differences between sites
    params.E_diff = zeros(lattice.m); %DEFAULT Site Energies
    for i = 1:lattice.m
        for j = 1:lattice.m

            % the energy difference between site i and site j
            params.E_diff(i,j) = params.E(j)-params.E(i);

        end
    end

    % Calculate the lambda matrix
    % lambda(1,j) is the ratio between site j to site i
    Gamma = @(e) exp(-e/(params.kB*params.T)); 
    params.lambda = Gamma(params.E_diff); % The ratio of site(energy state) concentrations given thermodynamically
    
end

% Calculate the relative site concentrations
% Eg. : for two sites, [c1 c2] = [c1 lambda*c1]
params.c = params.lambda(1,:)./sum(params.lambda(1,:)); %Normalized site concentration vectors

% Calculate the lower triangular part of tau matrix
% tau as input is an mxm upper triangular matrix.
% Theory: tau(j,i) = tau(i,j)*lambda(i,j)
for i = 1:lattice.m
    for j = (i+1):lattice.m
        params.tau(j,i) = params.tau(i,j)*params.lambda(i,j);
    end
end

%% Calculating each individual jump rate

params.jump_rates = zeros(lattice.m,lattice.m,max(max(lattice.singles)));
for i = 1:lattice.m
    for j = 1:lattice.m
        for k = 1:lattice.singles(i,j)
            params.jump_rates(i,j,k) = 1/(params.tau(i,j)*lattice.singles(i,j));
        end
    end
end

end
