function [out_mat_tn] = update_direction(dist,dist_xyz,in_mat_tn,in_mat_to)
% Input position and angle information for timepoint (t_old)
% into this function
%
% for consistency and atan2 and acos functions 
% theta is always between 0 and 2pi
% and phi is always between 0 and pi
%
% August 25, 2014

global L N r_c r_e r_a r_o X Y Z Theta eta alpha beta PHI

D_ae = 0.25/(r_a-r_e);

for cell = 1:N

    % cells with distance < radius
    C_e = dist(:,cell)<r_o;
    C_e_phi = dist(:,cell)<r_o & in_mat_to(:,PHI)~=0 & in_mat_to(:,PHI)~=pi;
    
    C_c  = dist(:,cell)<r_c;
    C_ca = (r_c<dist(:,cell)) & (dist(:,cell)<r_a);
    C_ao = (r_a<dist(:,cell)) & (dist(:,cell)<r_o);
    
    % Find average angle of neighbors' (and cell's own) migration
    N_nbors = nnz(C_e);   % Number of cells in r_e
    
    PHI_con = sum(cos(in_mat_to(C_e,PHI)));             % contribution of phi
    Theta_con = sum(exp(in_mat_to(C_e_phi,Theta)*1i));  % contribution of theta (excluding cells moving straight "up" (phi==0) or straight "down" (phi==pi))  
    
    %now find the argument of Z, only works if angles are between (-pi,pi)
    X_val1 = real(Theta_con);
    Y_val1 = imag(Theta_con);
    Z_val1 = PHI_con;
    
    % Normalize distance between cells and multiply based on Lennard-Jones(accounts for PBC)
    clear norm_c norm_ca norm_co;   % Contain distance in components (xyz) for cells within distance intervals (i.e. (0,r_c), (r_c,r_a), (r_a,r_o))
    
    C_c(cell) = 0;
    norm_c = (-10*L).*normr(dist_xyz(C_c,:,cell));
    
    if isempty(norm_c)          % if no cell is at <r_c, make it [0 0 0] to prevent NaN
        norm_c = zeros(1,3);
    end
    
    norm_ca(:,X:Z) = normr(dist_xyz(C_ca>0,X:Z,cell));    % D_ae = 1/4 * 1/(r_a - r_e)
    for dim = X:Z
        norm_ca(:,dim) = (((dist(C_ca,cell))-r_e)*D_ae).*norm_ca(:,dim);
    end
    
    if isempty(norm_ca)
        norm_ca = zeros(1,3);
    end
    
    norm_ao = normr(dist_xyz(C_ao,:,cell));
    
    if isempty(norm_ao)
        norm_ao = zeros(1,3);
    end
    
    % Add forces (attraction/repulsion) to affect angle
    X_val2 = sum(norm_c(:,X)) + sum(norm_ca(:,X)) + sum(norm_ao(:,X));
    Y_val2 = sum(norm_c(:,Y)) + sum(norm_ca(:,Y)) + sum(norm_ao(:,Y));
    Z_val2 = sum(norm_c(:,Z)) + sum(norm_ca(:,Z)) + sum(norm_ao(:,Z));
    
    % angles indicating direction of random unit vector u
    mag_noise = N_nbors*eta;
    ang_theta = 2*pi*rand; 
    ang_phi = acos(2*rand-1);
    
    % Add contributions of direction, forces, noise to new direction (theta & phi)
    X_val_angle = alpha*X_val1 + beta*X_val2 + mag_noise*cos(ang_theta)*sin(ang_phi);
    Y_val_angle = alpha*Y_val1 + beta*Y_val2 + mag_noise*sin(ang_theta)*sin(ang_phi);
    Z_val_angle = alpha*Z_val1 + beta*Z_val2 + mag_noise*cos(ang_phi);
    
    R_val_angle = realsqrt((X_val_angle^2 + Y_val_angle^2 + Z_val_angle^2));
    
    theta_avg = atan2(Y_val_angle,X_val_angle);
    %map to 0 to 2pi
    if theta_avg < 0
        theta_avg = theta_avg + 2*pi;
        
    elseif theta_avg == 0   % undefined case when Y_val_angle==0 and X_val_angle == 0
        theta_avg = in_mat_to(cell,Theta);     
    end
    
    phi_avg = acos(Z_val_angle/R_val_angle);
    if isnan(phi_avg) == 1 
        phi_avg = in_mat_to(cell,PHI); 
    end
    
    
    if isnan(theta_avg) == 1 
        theta_avg = in_mat_to(cell,Theta); 
        phi_avg = in_mat_to(cell,PHI);       
    end
    
    in_mat_tn(cell,Theta) = theta_avg;
    in_mat_tn(cell,PHI) = phi_avg;
end
  
out_mat_tn = in_mat_tn;

end


function B = normr(A)
    % Get the Euclidean norm (2-norm) of each row in matrix A
    rowNorms = sqrt(sum(A.^2, 2));
    
    % Avoid division by zero by setting zero norms to 1
    rowNorms(rowNorms == 0) = 1;
    
    % Normalize each row by dividing the row by its norm
    B = bsxfun(@rdivide, A, rowNorms);
end
