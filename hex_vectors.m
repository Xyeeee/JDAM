% Copyright (c) 2020, Xy Wang and Nadav Avidor.
% All rights reserved.
% This file is part of the JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.


%% Configure a hexagonal surface, Site order: T B_1 B_2 B_3 H_1 H_2

a1 = lattice.a*[1,0]; % The first lattice vector
a2 = lattice.a*[-1/2,sqrt(3)/2]; % The second lattice vector
z = [0,0]; % Zero Vector as a placeholder

%% Single jump vectors k
% The maximum number of distinct single jumps from site of type i to site of type
% j is 6, hence the third index k in L_ijk goes from 1 to 6, for those (i,j) pairs
% with less than 6 jumps, we pad the jump vectors with zero vectos z =
% [0,0]

% Default N matrix: Single Jumps only
lattice.singles = [0,2,2,2,3,3;...
                   0,0,2,2,1,1;...
                   0,0,0,2,1,1;...
                   0,0,0,0,1,1;...
                   0,0,0,0,0,3;...
                   0,0,0,0,0,0];
lattice.singles = lattice.singles + lattice.singles' + diag([6,0,0,0,0,0]);

lattice.s = zeros(6,6,6,2);

lattice.s(1,1,:,:) = [-a1;a1;-a2;a2;a1+a2;-a1-a2];%6
lattice.s(1,2,:,:) = [a1/2;-a1/2;z;z;z;z];%2
lattice.s(1,3,:,:) = [a2/2;-a2/2;z;z;z;z];%2
lattice.s(1,4,:,:) = [(a1+a2)/2;-(a1+a2)/2;z;z;z;z];%2
lattice.s(1,5,:,:) = [-a1/2+1/3*(a1/2+a2);a1/2+1/3*(a1/2+a2);-2/3*(a1/2+a2);z;z;z];%3
lattice.s(1,6,:,:) = [-a1/2-1/3*(a1/2+a2);a1/2-1/3*(a1/2+a2);2/3*(a1/2+a2);z;z;z];%3
lattice.s(2,2,:,:) = [z;z;z;z;z;z];%2
lattice.s(2,3,:,:) = [a2/2;-a2/2;z;z;z;z];%2
lattice.s(2,4,:,:) = [(a1+a2)/2;-(a1+a2)/2;z;z;z;z];%2
lattice.s(2,5,:,:) = [1/3*(a1/2+a2);z;z;z;z;z];%1
lattice.s(2,6,:,:) = [-1/3*(a1/2+a2);z;z;z;z;z];%1
lattice.s(3,3,:,:) = [z;z;z;z;z;z];%2
lattice.s(3,4,:,:) = [-a1/2;a1/2;z;z;z;z];%2
lattice.s(3,5,:,:) = [-a1/2+1/3*(a1/2+a2)-a2/2;z;z;z;z;z];%1
lattice.s(3,6,:,:) = [2/3*(a1/2+a2)-a2/2;z;z;z;z;z];%1
lattice.s(4,4,:,:) = [z;z;z;z;z;z];%2
lattice.s(4,5,:,:) = [a1/2+1/3*(a1/2+a2)-(a1+a2)/2;z;z;z;z;z];%1
lattice.s(4,6,:,:) = [2/3*(a1/2+a2)-(a1+a2)/2;z;z;z;z;z];%1
lattice.s(5,6,:,:) = [-a1/2+1/3*(a1/2+a2);a1/2+1/3*(a1/2+a2);-2/3*(a1/2+a2);z;z;z];%3

%% Complete Lattice Jump Vectors

lattice.m = size(params.tau,1);

% The jump vector from i to j is opposite in direction to that from j to i
% and same in magnitude
for i = 1:6
    for j = i+1:6
        lattice.s(j,i,:,:) = lattice.s(i,j,:,:).*(-1);
    end
end


