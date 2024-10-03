function cell_conc = initialize_rcn()
% holds parameters for RCN (as global variables), and initializes molecular
% concentrations
%
% January 29, 2015

global N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RCN Parameters
% found in m-files: initialize_rcn.m, solve_rcn.m, rcn_functions
% (within solve_rcn)

WNT0 = 28.062;      APC0 = 100;         TCF0 = 15;          AxGSK0 = 0.02;
K1 = 6;             K2 = 1/12;          k4 = 0.26;          k5 = 0.13;
k3 = 0.091;         km3 = 0.91;         K6 = 100;           k7 = 210;
K8 = 1000;          v9 = 0.599976;      k10 = 0.00026;      K11 = 30;
Tmax = 0.005946;    Krna = 10;          yps = 3;            Pmax = 0.025;
k13 = 0.11781;      v15 = 0.1;          k16 = 0.02;         k17 = 0.1;
k18 = 0.015;        k19 = 0.0924;       v20 = 0.000024;     k21 = 0.00006412;
k22 = 0.00231049;   k23 = 0.006001;     k25 = 0.00025;      VtG = 1.1;
Km = 0.00886;       kaj = 0.0001068;    kdaj = 0.0002136;

var = 22;                       % Number of variables in RCN (molecules)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_conc = zeros(N,var,2);     % Keep concentrations of only the last two times steps for every cell

% Steady-states for WNT0 = 28.062 while  respecting all conservation equations:
X2o = 1.29122;
X3o = 0.014274;
X1o = WNT0 - X2o - X3o;

X4o = K2*X3o/X2o;
X5o = 0.00184245;

X6o = 0.000106402;
X7o = 86.5796;
X8o = 13.4204;

X9o = 0.0028559;
X10o = 155.006;
X11o = 2.43235;
X12o = 12.5676;

X13o = 0.0335629;
X14o = 0.00908087;
X15o = 4.78284;
X16o = K1*X2o/X1o;

X17o = 0.374298;
X18o = 0.259727;
X19o = 0.096;

SigERo = 0.443436;
SigMo = 0.443436;
SigERCo = 0.443436;

% Define initial concentrations
cell_conc(:,:,2) = repmat([X1o X2o X3o X4o X5o X6o X7o X8o X9o X10o X11o X12o ...
        X13o X14o X15o X16o X17o X18o X19o SigERo SigMo SigERCo],N,1);
