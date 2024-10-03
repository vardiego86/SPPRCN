function [cell_mat] = initialize()
%for consistency and atan2 and acos funcitons 
%theta is always between 0 and 2pi
% and phi is always between 0 and pi
%
% August 1, 2014

global L N T X Z t_p Theta PHI ID

cell_mat = zeros(N,7,T);

% Position
cell_mat(:,X:Z,1) = L*rand(N,3);
    
% Direction
cell_mat(:,Theta,1) = 2*pi*rand(N,1);   % Theta- between 0 and 2*pi
cell_mat(:,PHI,1) = pi*rand(N,1);       % PHI -between  0 and pi

% time point
cell_mat(:,t_p,:) = repmat(1:1:T,N,1);

% cell ID
cell_mat(:,ID,:) = repmat((1:1:N)',1,T); 

