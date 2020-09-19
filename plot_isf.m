% Copyright (c) 2020, Yaqing Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

function plot_isf(Azmth,azmthStr,dK)
%% Loading surface configurations
load('surface');

if isempty(azmthStr), azmthStr = cell(zeros(size(Azmth,1))); end

% Preparing to plot
mrkrs = 'o+ds>*'; % Data sequence markers
figure; % Creating a figure

% Computing the A matrix to extract its eigenvalues(decay constants) and
% eigenvectors(later processed to be intensities)
for x = 1: size(Azmth,1)
    
    tSE = 0:1:300;
    [Wp_all,D_abs,ISF] = calc_ISF(Azmth(x,:),dK,lattice,jump_rates,'tSE',tSE);
    
    if size(lattice.tau,1) > 1 && ~isinf(lattice.tau(1,2))
        scaleTau = lattice.tau(2,1)^(-1);
        ylabelScale = '/$\tau_{21}^{-1}$';
        titleScale = ' SCALED';
    else
        scaleTau = 1;
        ylabelScale = '';
        titleScale = '';
    end
    
    D_norm = D_abs/scaleTau;
    
    % Test mode: zero all decay rates which correspond to near zero
    % intensity
    if 1 && max(D_norm(:))-min(D_norm(:)) > 100
        for i = 1:lattice.m
            nearZero = abs(Wp_all)<1e-4;
            Wp_all(nearZero) = 0;
            D_norm(nearZero) = 0;
        end
    end
    
    
    % Plotting 
    subplot(size(Azmth,1),3,(x-1)*3+1)
    xlabel('dK/(2pi/a)','Interpreter','Latex','FontSize',14)
    ylabel(['Decay' ylabelScale],'Interpreter','Latex','FontSize',14)    
    title([azmthStr{x} titleScale ' Decay constants'],'Interpreter','Latex','FontSize',14)
    hold on
    for i = 1:lattice.m
        plot(dK,D_norm(i,:),mrkrs(i))        
    end
    hold off
    
    subplot(size(Azmth,1),3,(x-1)*3+2)
    xlabel('dK/(2pi/a)','Interpreter','Latex','FontSize',14)
    ylabel('Intensity','Interpreter','Latex','FontSize',14)
    title([azmthStr{x} 'Intensities'],'Interpreter','Latex','FontSize',14)        
    hold on
    for i = 1:lattice.m
        plot(dK,Wp_all(i,:),mrkrs(i))        
    end
    hold off
    
    subplot(size(Azmth,1),3,(x-1)*3+3)
    xlabel('tSE [ps]','Interpreter','Latex','FontSize',14)
    ylabel('ISF','Interpreter','Latex','FontSize',14)
    title([azmthStr{x} 'Intensities'],'Interpreter','Latex','FontSize',14)        
    hold on
    for i=1:10:length(dK)
        plot(tSE,ISF(i,:));
    end
    hold off
    
end

end % end plot_isf



