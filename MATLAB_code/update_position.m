function [out_mat_tn] = update_position(in_mat_tn,in_mat_to,junc)
% Input position and angle information for timepoint (t_old)
% into this function
%
% scales speed based on adherens junctions
%
% August 6, 2014

global L X Y Z Theta PHI speed

% Scaling factor based on total AJs formed by cell

drop = 120000;      % Affects sharpness of drop in speed (in scaling factor from 1 to 0)
range = 4;          % Shifts drop on range of [AJ]

junc(junc==0) = nan;
aj_scaling = (-1./(1+exp(-((nanmean(junc,2))*drop-range))))+1;   % Sigmoidal function relating number of junctions to scaling of speed ( y(x) = -1/(1 + e-x) + 1 )
aj_scaling(isnan(aj_scaling)) = 1;

aj_scaling = aj_scaling.*speed;

% New position (XYZ) based on NEW theta
in_mat_tn(:,X) = in_mat_to(:,X) + aj_scaling.*cos(in_mat_tn(:,Theta)).*sin(in_mat_tn(:,PHI));
in_mat_tn(:,Y) = in_mat_to(:,Y) + aj_scaling.*sin(in_mat_tn(:,Theta)).*sin(in_mat_tn(:,PHI));
in_mat_tn(:,Z) = in_mat_to(:,Z) + aj_scaling.*cos(in_mat_tn(:,PHI));


% Accounting for PBC (periodic boundary conditions)
% X
Xf = in_mat_tn(:,X) > L;
in_mat_tn(Xf,X) = in_mat_tn(Xf,X) - L;

Xf = in_mat_tn(:,X) < 0;
in_mat_tn(Xf,X) = in_mat_tn(Xf,X) + L;

% Y
Yf = in_mat_tn(:,Y) > L;
in_mat_tn(Yf,Y) = in_mat_tn(Yf,Y) - L; 

Yf = in_mat_tn(:,Y) < 0;
in_mat_tn(Yf,Y) = in_mat_tn(Yf,Y) + L;  

% Z
Zf = in_mat_tn(:,Z) > L;
in_mat_tn(Zf,Z) = in_mat_tn(Zf,Z) - L; 

Zf = in_mat_tn(:,Z) < 0;
in_mat_tn(Zf,Z) = in_mat_tn(Zf,Z) + L;  


out_mat_tn = in_mat_tn;