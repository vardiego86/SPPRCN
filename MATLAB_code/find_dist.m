function [dist,dist_xyz] = find_dist(cell_mat_t)
% find distance between every pair of cells at a particular time
%
% August 1, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global L N X Y Z

dist = zeros(N,N);          % dist: array to calculate distance from each other cell   
dist_xyz = zeros(N,3,N);      % xyz distance for every cell

for cell = 1:N
   
    dist_xyz(:,X,cell) = cell_mat_t(:,X)-cell_mat_t(cell,X);
    dist_xyz(:,Y,cell) = cell_mat_t(:,Y)-cell_mat_t(cell,Y);
    dist_xyz(:,Z,cell) = cell_mat_t(:,Z)-cell_mat_t(cell,Z);

    % Accounting for PBC (periodic boundary conditions)
    %All greater than L/2
    Xf = dist_xyz(:,X,cell) > L/2;
    dist_xyz(Xf,X,cell) = dist_xyz(Xf,X,cell) - L;
    Yf = dist_xyz(:,Y,cell) > L/2;
    dist_xyz(Yf,Y,cell) = dist_xyz(Yf,Y,cell) - L;
    Zf = dist_xyz(:,Z,cell) > L/2;
    dist_xyz(Zf,Z,cell) = dist_xyz(Zf,Z,cell) - L;
    
    %All less than L/2
    Xf = dist_xyz(:,X,cell) <= -L/2;
    dist_xyz(Xf,X,cell) = dist_xyz(Xf,X,cell) + L;
    Yf = dist_xyz(:,Y,cell) <= -L/2;
    dist_xyz(Yf,Y,cell) = dist_xyz(Yf,Y,cell) + L;
    Zf = dist_xyz(:,Z,cell) <=-L/2;
    dist_xyz(Zf,Z,cell) = dist_xyz(Zf,Z,cell) + L;
    
    dist(cell,:) = realsqrt(sum(dist_xyz(:,X:Z,cell).^2,2));

end