% Copyright (c) 2020, Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

function [Wp_all,D_abs,ISF] = plot_isf(lattice,params)

global scaleDecayInPlot

if isempty(params.azmthStr), params.azmthStr = cell(zeros(size(params.Azmth,1))); end

% Preparing to plot
mrkrs = 'o+ds>*'; % Data sequence markers
figure; % Creating a figure

% Computing the A matrix to extract its eigenvalues(decay constants) and
% eigenvectors(later processed to be intensities)
for x = 1: size(params.Azmth,1)
    
    [Wp_all,D_abs,ISF] = calc_ISF(params.Azmth(x,:),params,lattice);
    
    
    % ==== Scale and Plot ==== %
    
    if scaleDecayInPlot && size(params.tau,1) > 1 && ~isinf(params.tau(1,2))
        scaleTau = params.tau(2,1)^(-1);
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
    
    % Scale dK if requested
    if params.dKscale == 2*pi && lattice.a == 1
        dKscale = params.dKscale;
        dKStr = '/(2pi/a)';
    else
        dKStr = '';
        dKscale = 1;
    end
    
    % Plotting
    nPlots = 3;
    
    subplot(size(params.Azmth,1),nPlots,(x-1)*nPlots+1)
    xlabel(['dK' dKStr],'Interpreter','Latex','FontSize',14)
    ylabel(['Decay' ylabelScale],'Interpreter','Latex','FontSize',14)    
    title([params.azmthStr{x} titleScale ' Decay constants'],'Interpreter','Latex','FontSize',14)
    hold on
    for i = 1:lattice.m
        plot(params.dK/dKscale,D_norm(i,:),mrkrs(i))        
    end
    hold off
    
    subplot(size(params.Azmth,1),nPlots,(x-1)*nPlots+2)
    xlabel(['dK' dKStr],'Interpreter','Latex','FontSize',14)
    ylabel('Intensity','Interpreter','Latex','FontSize',14)
    title([params.azmthStr{x} 'Intensities'],'Interpreter','Latex','FontSize',14)        
    hold on
    for i = 1:lattice.m
        plot(params.dK/dKscale,Wp_all(i,:),mrkrs(i))        
    end
    hold off
    
    subplot(size(params.Azmth,1),nPlots,(x-1)*nPlots+3)
    xlabel('tSE [ps]','Interpreter','Latex','FontSize',14)
    ylabel('ISF','Interpreter','Latex','FontSize',14)
    title([params.azmthStr{x} 'Intensities'],'Interpreter','Latex','FontSize',14)        
    hold on
    for i=1:1:length(params.dK)
        plot(params.tSE,ISF(i,:));
    end
    hold off
    
end

end % end plot_isf



