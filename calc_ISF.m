% Copyright (c) 2020,Yaqing Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

function [Wp_all,D_abs,ISF] = calc_ISF(K_base,dK,lattice,jump_rates,varargin)

% parse varargin
prsdArgs = inputParser;   % Create instance of inputParser class.            
prsdArgs.addParameter('tSE',0, @isnumeric); % in ps
prsdArgs.parse(varargin{:});            
tSE = prsdArgs.Results.tSE;            

D_abs = zeros(lattice.m,length(dK)); % decay constant vector to be computed
Wp_all = zeros(lattice.m,length(dK)); % intensity(normalized) vector to be computed
ISF = zeros(length(dK),length(tSE));

% ====== Calculate Dephasing Rates and Intensities ======
for q = 1:length(dK)
    
    %% Calculate the A matrix for current dK vector K_base.*dK(q)
    A = zeros(lattice.m);
    for i = 1:lattice.m
        for j = 1:lattice.m
            rum = 0;
            for first = 1:lattice.singles(i,j)
                jmpVec = squeeze( lattice.s(i,j,first,:) );
                rum = rum + jump_rates(j,i,first)*exp(-1i*dot(K_base.*dK(q),jmpVec)*2*pi);
            end
            if i==j
                for tmpj = 1:lattice.m
                    % Take into account intra lattice jumps, unless its not
                    % defined as part of the lattice jump vectors
                    if i == tmpj && lattice.singles(i,tmpj) == 0
                        continue
                    end
                    rum = rum - 1/lattice.tau(i,tmpj);
                end
            end
            A(i,j) = rum;
        end
    end
    
    % Transfrom A into a Hermitian Matrix with T
    T = diag(sqrt(1./lattice.c));

    % Get B the Hermitian Matrix
    B = T*A/T;


    % Solve for eigenvalues and eigenvectors of A
    [V,D] = eig(B);

    D = diag(D);
    D_abs(:,q) = abs(D);

    for p = 1:lattice.m
        Wp = 0;
        for i = 1:lattice.m
            Wp = Wp + V(i,p)*sqrt(lattice.c(i));
        end
        Wp_all(p,q) = abs(Wp)^2;
        
        ISF(q,:) = ISF(q,:) + Wp_all(p,q).*exp(-D_abs(p,q).*tSE);
    end
   
end

end % end calc_ISF