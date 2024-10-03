%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 3D multiscale model: SPP (based of 2D Gregoire et al. Moving and staying 
% together witout a leader. Physica D 181: 157-170 (2003)) + RCN - reduced
% model
%
% Series version (solves RCN for groups of cells in series instead of 
% parallel).
%
% Saves concentrations for all time steps, as well
% as average [AJ] for each cell and corresponding sigma values.
% 
% Yasha Sharma and Diego A. Vargas
% Boston University
% March 12, 2015
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

%%%%%%%%%%%%%%%

alphai = 3;
betai = 3;

display(alphai);
display(betai);

global L N T r_c r_e r_a r_o X Y Z t_p Theta ID eta alpha beta speed PHI a_neigh

% 0) Initialize
%identifiers for columns of main matrix (cell_mat)
X =     1;
Y =     2;
Z =     3; 
t_p=    4;
Theta = 5; 
PHI =   6;
ID =    7;

T = 5000;           % Duration of simulation (No. steps, SPP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETERS

L = 1280;           % Dimension world
N = 512;            % Number of particles

% Radii of interaction
r_c = 8;            % core repulsion radius
r_e = 20;           % equilibrium (preferred) radius
r_a = 32;           % limit of attraction radius
r_o = 36;           % limit of interaction (changed from paper to account for lack of Voronoi tesselation)

a_neigh = 12;            % allowed max number of neighbors per cell

% Noise parameter
eta = 1;

% Relative importance of forces
alpha = alphai;
beta = betai;

speed = 2;       % spatial units/t_step

% RCN number of variables (var) and steps (tr) values are established in
% solve_rcn.m to avoid using global variables in a parallel environment.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mol_conc = zeros(N,22,T);       % Save concentrations for last 5000 time points,
aj_conc = zeros(N,T);           % including average aj concentration per cell,
sigAJ_conc = zeros(N,T);        % and adhesivity factor of AJs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Assign initial position and direction
cell_mat = initialize();

[cell_dist,cell_dist_xyz] = find_dist(cell_mat(:,:,1));

% 2) Initialize RCN:
cell_conc = initialize_rcn();                             % initial concentrations
cell_junc = initialize_junc(cell_dist,cell_conc(:,:,2));    % calculate AJs based on neighbor relations

SigAJo = 0.443436;
sigma_aj = SigAJo.*ones(N,1);

%allocate memory for cell_conc in parfor calculations
cell_conc_old = zeros(22,N);
cell_conc_new = zeros(22,N);

%3) Update direction (Theta,Phi) and position (XYZ)
for t = 2:T

    if ~mod(t,5)
        display(t);
    end

    t_new = t;
    t_old = t_new - 1;

    % a) SPP alone: Update direction (purely based on neighbor relations)
    cell_mat(:,:,t_new) = update_direction(cell_dist,cell_dist_xyz,cell_mat(:,:,t_new),cell_mat(:,:,t_old));

    % b) RCN alone: Solve RCN model
    cell_conc(:,:,1) = cell_conc(:,:,2);
    
    % transpose to pass to solve_rcn
    cell_conc_old = cell_conc(:,:,1)';
   
    for celli = 1:N
        cell_conc_new(:,celli) = solve_rcn_ser(cell_conc_old(:,celli));
    end
    cell_conc(:,:,2) = cell_conc_new';
    
    mol_conc(:,:,t) = cell_conc(:,:,2);

    % c) SPP affecting RCN: Update [AJ] and (E-cad/B-cat)M
    [cell_dist,cell_dist_xyz] = find_dist(cell_mat(:,:,t_old));

    cell_junc(:,:,1) = cell_junc(:,:,2);
    [cell_conc(:,:,2),cell_junc] = find_junc(cell_dist,cell_conc(:,:,2),cell_junc,sigma_aj);    % find AJs and update concentrations
    
    aj_conc(:,t) = nanmean(cell_junc(:,:,2),2);                                 % save average AJ concentrations

    [sigma_aj,cell_conc(:,:,2)] = update_sigma(sigma_aj,cell_conc(:,:,2),cell_junc(:,:,2));
    
    sigAJ_conc(:,t) = sigma_aj;
    
    % d) RCN affecting SPP: Update positions based on direction and [AJ](scaling speed)
    cell_mat(:,:,t_new) = update_position(cell_mat(:,:,t_new),cell_mat(:,:,t_old),cell_junc(:,:,2));
    
    % Save statistics in .mat file every 1000 steps
    if ~mod(t,1000)
        filename = strcat('Alpha_',num2str(alphai),'_Beta_',num2str(betai),'_step_',num2str(t),'.mat');
        save(filename,'cell_mat','mol_conc','aj_conc','sigAJ_conc','-v7.3');
    end

end


