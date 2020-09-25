% Copyright (c) 2020, Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

function [Wp_all,D_abs,ISF] = calc_ISF(K_base,params,lattice)           

if params.Psurvival ~= -1
    % Assume an explicit model for MULTIPLE jumps on hollow sites
    % (Implemented as an interim solution until multiple jumps are
    % implemented for the general case)
    if isfield(lattice,'m') && lattice.m == 2
    else
        lattice.m = 2;
        params.c = [1 params.lambda]/(1+params.lambda);
    end
end

D_abs = zeros(lattice.m,length(params.dK)); % decay constant vector to be computed
Wp_all = zeros(lattice.m,length(params.dK)); % intensity(normalized) vector to be computed
ISF = zeros(length(params.dK),length(params.tSE));

% ====== Calculate Dephasing Rates and Intensities ======
for q = 1:length(params.dK)
    
    %% Calculate the A matrix for current params.dK vector K_base.*params.dK(q)
    
    % Calculate the general case, unless requested otherwise
    % (the general case is currently indicated by lattice.Psurvival == -1)
    if params.Psurvival == -1
        
        A = zeros(lattice.m);
        for i = 1:lattice.m
            for j = 1:lattice.m
                rum = 0;
                for first = 1:lattice.singles(i,j)
                    jmpVec = squeeze( lattice.s(i,j,first,:) );
                    rum = rum + params.jump_rates(j,i,first)*exp(-1i*dot(K_base.*params.dK(q),jmpVec));
                end
                if i==j
                    for tmpj = 1:lattice.m
                        % Take into account intra lattice jumps, unless its not
                        % defined as part of the lattice jump vectors
                        if i == tmpj && lattice.singles(i,tmpj) == 0
                            continue
                        end
                        rum = rum - 1/params.tau(i,tmpj);
                    end
                end
                A(i,j) = rum;
            end
        end
        
    else
        % if lattice.Psurvival is not -1, calc A from multiple jumps up to
        % 4th order. This assumes that the user have set all the relevant
        % parameters in the 'lattice' structure: tau, lambda, Psurvival
        
        % if tau and lambda are 2x2 matrix (which will be the case if the 
        % package was used to generate it out of the 'general' case), reduce it.
        if numel(params.tau)>1, params.tau = params.tau(1,2); end
        if numel(params.lambda)>1, params.lambda = params.lambda(1,2); end
        if numel(params.c) ~= 2, error('C must be of size 2'); end
        
        A = multiple_hollow_A(lattice.a,params.Psurvival,params.tau,params.lambda,K_base,params.dK(q));
        
    end
    
    % Transfrom A into a Hermitian Matrix with T
    T = diag(sqrt(1./params.c));

    % Get B the Hermitian Matrix
    B = T*A/T;


    % Solve for eigenvalues and eigenvectors of A
    [V,D] = eig(B);

    D = diag(D);
    D_abs(:,q) = abs(D);

    for p = 1:lattice.m
        Wp = 0;
        for i = 1:lattice.m
            Wp = Wp + V(i,p)*sqrt(params.c(i));
        end
        Wp_all(p,q) = abs(Wp)^2;
        
        ISF(q,:) = ISF(q,:) + Wp_all(p,q).*exp(-D_abs(p,q).*params.tSE);
    end
   
end

end % end calc_ISF