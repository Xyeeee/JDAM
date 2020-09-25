% Copyright (c) 2020, Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

function A = multiple_hollow_A(a,s,tau,lambda,K_base,dK)
%MULTIPLE_HOLLOW_A calculates the 'A' matrix for jump diffusion on hollow
%sites. It considers currently jumps to 10th nearest neighbour, and assumes
%a survival probability model in the form of a geometric series.
%
% Inputs:
%        a - lattice constant
%        s - survival probability for staing in the mobile state and not adsorbing in the nearest neighbour
%        tau - residence time at the higher energy site, [ps]
%        lambda - a measure of the energetic inequivalence between fcc and hcp sites
%        azimuthal direction (normally [0 1] or [1 0], but can be anything less trivial
%        dK - momentum transfer in Angstrom
%


%% User parameters
% s = 1;
% tau = 3; lambda = 3;
% K_base = [1,0];
% dK = 1;
% a = 2.71; % Ru0001

dKvec = K_base*dK;
p = [1 s s^2 s^3 s^4 s^5 s^6 s^7 s^8 s^9]; p = p/sum(p);

a1 = a*[1,0]; % The first lattice vector
a2 = a*[-1/2,sqrt(3)/2]; % The second lattice vector

%% single jumps
v1_12 = [-a1/2+1/3*(a1/2+a2);a1/2+1/3*(a1/2+a2);-2/3*(a1/2+a2)];

%% Double jumps
v2_11 = [];
for i=1:size(v1_12,1)
    v2_11 = [v2_11;v1_12(i,:)-v1_12];
end
v2_11 = filter_jump_by_length(v2_11,v1_12);

%% Triple jumps
v3_12 = [];
for i=1:size(v2_11,1)
    v3_12 = [v3_12;v2_11(i,:)+v1_12];
end
v3_12 = filter_jump_by_length(v3_12,v2_11);

%% Quatruple jumps
v4_11 = [];
for i=1:size(v3_12,1)
    v4_11 = [v4_11;v3_12(i,:)-v1_12];
end
v4_11 = filter_jump_by_length(v4_11,v3_12);

%% Quituple jumps
v5_12 = [];
for i=1:size(v4_11,1)
    v5_12 = [v5_12;v4_11(i,:)+v1_12];
end
v5_12 = filter_jump_by_length(v5_12,v4_11);

%% Sextuple jumps
v6_11 = [];
for i=1:size(v5_12,1)
    v6_11 = [v6_11;v5_12(i,:)-v1_12];
end
v6_11 = filter_jump_by_length(v6_11,v5_12);

%% Septuple jumps
v7_12 = [];
for i=1:size(v6_11,1)
    v7_12 = [v7_12;v6_11(i,:)+v1_12];
end
v7_12 = filter_jump_by_length(v7_12,v6_11);

%% Octuple jumps
v8_11 = [];
for i=1:size(v7_12,1)
    v8_11 = [v8_11;v7_12(i,:)-v1_12];
end
v8_11 = filter_jump_by_length(v8_11,v7_12);

%% Nonuple jumps
v9_12 = [];
for i=1:size(v8_11,1)
    v9_12 = [v9_12;v8_11(i,:)+v1_12];
end
v9_12 = filter_jump_by_length(v9_12,v8_11);

%% Decuple jumps
v10_11 = [];
for i=1:size(v9_12,1)
    v10_11 = [v10_11;v9_12(i,:)-v1_12];
end
v10_11 = filter_jump_by_length(v10_11,v9_12);

%%

% %% Plot
% figure
% plot(v1_12(:,1),v1_12(:,2),'ko',v2_11(:,1),v2_11(:,2),'ro',...
%      v3_12(:,1),v3_12(:,2),'bo',v4_11(:,1),v4_11(:,2),'mo');
 
A = calc_A(p,tau,lambda,dKvec,v1_12,v2_11,v3_12,v4_11,v5_12,v6_11,v7_12,v8_11,v9_12,v10_11);


%% Auxilary Functions

function [v,r] = filter_jump_by_length(v,v_origin)

% r = [r_min r_max], the range within to accept lengths of vectors
r_min = min(sqrt(sum(v_origin.^2,2)));
r_max = max(sqrt(sum(v_origin.^2,2)));
r = [r_min 2*r_max];

l = sqrt(sum(v.^2,2));

ind = r(1)<l & l<r(2);
v=v(ind,:);
l = l(ind);

% Take into account a situation were two paths lead to the same site
l_unique = unique(round(l,5));
if length(l_unique)==2
    
end

v_new = [];

end

function A = calc_A(p,tau,lambda,dKvec,v1_12,v2_11,v3_12,v4_11,v5_12,v6_11,v7_12,v8_11,v9_12,v10_11)
 % Calc A
 
A_11_2th    = sum(p(2)/tau*exp(-1i*dKvec*v2_11'))/size(v2_11,1);
A_11_4th    = sum(p(4)/tau*exp(-1i*dKvec*v4_11'))/size(v4_11,1);
A_11_6th    = sum(p(6)/tau*exp(-1i*dKvec*v6_11'))/size(v6_11,1);
A_11_8th    = sum(p(8)/tau*exp(-1i*dKvec*v8_11'))/size(v8_11,1);
A_11_10th   = sum(p(10)/tau*exp(-1i*dKvec*v10_11'))/size(v10_11,1);

A_11 = A_11_2th + A_11_4th + A_11_6th + A_11_8th + A_11_10th - 1/tau;

A_22 = A_11/lambda;

A_12_1th   = sum(p(1)/tau*exp(-1i*dKvec*v1_12'))/size(v1_12,1);
A_12_3th   = sum(p(3)/tau*exp(-1i*dKvec*v3_12'))/size(v3_12,1);
A_12_5th = sum(p(5)/tau*exp(-1i*dKvec*v5_12'))/size(v5_12,1);
A_12_7th = sum(p(7)/tau*exp(-1i*dKvec*v7_12'))/size(v7_12,1);
A_12_9th = sum(p(9)/tau*exp(-1i*dKvec*v9_12'))/size(v9_12,1);

A_12 =     (A_12_1th + A_12_3th + A_12_5th + A_12_7th + A_12_9th)/lambda;

A_21 = conj(A_12_1th + A_12_3th + A_12_5th + A_12_7th + A_12_9th);

A = [A_11 A_12; A_21 A_22];
end



end