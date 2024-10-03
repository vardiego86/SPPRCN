function conc_out = solve_rcn_ser(conc_in)

% Diego A. Vargas
% Boston University
% February 5, 2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Differential equation solver for proposed model (Wnt-canonical pathway
% centered on N-glycosylation) based on Kirschner model(Le et al., PLoS Bio
% 2003)

% Uses MATLAB ode15i to solve DAEs (differential algebraic equations)

% NO ADHERENS JUNCTIONS IN SYSTEM

% Checks to see if negative values of concentrations are encountered after
% solving

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conc_in = conc_in';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RCN Parameters
% found in m-files: initialize_rcn.m, solve_rcn.m, rcn_functions
% (withing solve_rcn)

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

% time steps
tr = 10;        % How many rcn time step per spp time step

tstart = 0;
tend = tr-1;
n = tr;
tspan = linspace(tstart,tend,n);

X1 = 1;     X2 = 2;     X3 = 3;     X4 = 4;     X5 = 5;     
X6 = 6;     X7 = 7;     X8 = 8;     X9 = 9;     X10 = 10;
X11 = 11;   X12 = 12;   X13 = 13;   X14 = 14;   X15 = 15;   
X16 = 16;   X17 = 17;   X18 = 18;   X19 = 19;  
SigER = 20;     SigM = 21;      SigERC = 22;

% C initial conditions
X1o = conc_in(1,1);     X2o = conc_in(1,2);     X3o = conc_in(1,3);
X4o = conc_in(1,4);     X5o = conc_in(1,5);     X6o = conc_in(1,6);
X7o = conc_in(1,7);     X8o = conc_in(1,8);     X9o = conc_in(1,9);
X10o = conc_in(1,10);

X11o = conc_in(1,11);   X12o = conc_in(1,12);   X13o = conc_in(1,13);
X14o = conc_in(1,14);   X15o = conc_in(1,15);   X16o = conc_in(1,16);
X17o = conc_in(1,17);   X18o = conc_in(1,18);   X19o = conc_in(1,19);

SigERo = conc_in(1,20);
SigMo = conc_in(1,21);
SigERCo = conc_in(1,22);

% Fixing those for independent variables:
C0 = zeros(1,var);
fixed_C0 = zeros(1,var);

C0(X1) = X1o;           fixed_C0(X1) = 0;
C0(X2) = X2o;           fixed_C0(X2) = 0;
C0(X3) = X3o;           fixed_C0(X3) = 0;
C0(X4) = X4o;           fixed_C0(X4) = 1;
C0(X5) = X5o;           fixed_C0(X5) = 1;
C0(X6) = X6o;           fixed_C0(X6) = 1;
C0(X7) = X7o;           fixed_C0(X7) = 0;
C0(X8) = X8o;           fixed_C0(X8) = 0;
C0(X9) = X9o;           fixed_C0(X9) = 0;
C0(X10) = X10o;         fixed_C0(X10) = 1;
C0(X11) = X11o;         fixed_C0(X11) = 0;
C0(X12) = X12o;         fixed_C0(X12) = 0;
C0(X13) = X13o;         fixed_C0(X13) = 1;
C0(X14) = X14o;         fixed_C0(X14) = 1;
C0(X15) = X15o;         fixed_C0(X15) = 0;
C0(X16) = X16o;         fixed_C0(X16) = 0;
C0(X17) = X17o;         fixed_C0(X17) = 1;
C0(X18) = X18o;         fixed_C0(X18) = 1;
C0(X19) = X19o;         fixed_C0(X19) = 1;
C0(SigER) = SigERo;     fixed_C0(SigER) = 1;
C0(SigM) = SigMo;       fixed_C0(SigM) = 1;
C0(SigERC) = SigERCo;   fixed_C0(SigERC) = 1;

% dCdt initial conditions (guesses for dedic function)
dCdt0 = zeros(1,var);

dCdt0(X6) = km3*X4o - k3*X6o*X7o ;

dCdt0(X13) = Tmax*(X12o^yps/(Krna^yps + X12o^yps)) - k13*X13o;
dCdt0(X14) = Pmax*X13o - k19*X14o;

dCdt0(X15) = v15 - k17*X14o*X15o - k16*X15o;

dCdt0(X17) = v20 - k21*X17o;
dCdt0(X18) = k21*X17o - k22*X18o + k23*X19o;
dCdt0(X19) = k22*X18o - k23*X19o - k25*X19o;

dCdt0(SigER) = (v20/X17o)*((1 - (VtG*X14o/(Km+X14o))) - SigERo);
dCdt0(SigM) = (k21*X17o/X18o)*(SigERo - SigMo) + (k23*X19o/X18o)*(SigERCo - SigMo);
dCdt0(SigERC) = (k22*X18o/X19o)*(SigMo - SigERCo);

% System of equations to solve for dCdt0(X5) and dCdt0(X10):
% (1) -dX5/dt + dX8/dt + dX10/dt + dX12/dt = ...
% (2) dX5/dt + dX9/dt = ...

M = zeros(2,2);
M(1,1) = -1 ;
M(1,2) =  APC0*K8/(K8+X10o)^2 + TCF0*K11/(K11+X10o)^2 + 1 ;

M(2,1) = X10o/K6 + 1 ;
M(2,2) = X5o/K6 ;

c = zeros(2,1);
c(1,1) = -k4*X4o + k5*X5o - k7*X9o + v9 + k25*X19o - k10 ;
c(2,1) = k4*X4o + k5*X5o ;

y = M\c;

dCdt0(X5) = y(1);
dCdt0(X10) = y(2);

% some dependent variables:

dCdt0(X7) = -APC0*K8/(K8+X10o)^2*dCdt0(X10) ;
dCdt0(X8) = -dCdt0(X7) ;
dCdt0(X9) = (X10o/K6)*dCdt0(X5) + (X5o/K6)*dCdt0(X10) ;

dCdt0(X11) = TCF0*K11/(K11+X10o)^2*dCdt0(X10) ;
dCdt0(X12) = -dCdt0(X11) ;

% System of equations to solve for dCdt0(X3), dCdt0(X4) and dCdt0(X16):
% (1) -dX1/dt + dX3dt + dX4/dt + dX16 = ...
% (2) AxGSK0 = X3 + X4 + X5 + X6 + X9
% (3) dX3/dt(dX4/dt,dX16dt) = ..

O = zeros(3,3);
O(1,1) = 0 ;
O(1,2) = ((K1^2*K2^2+2*K1*K2*X16o*(K2+WNT0+X4o)+X16o^2*(K2^2+X4o^2+K2*(WNT0+2*X4o)))/(K1*K2+X16o*(K2+X4o))^2) ;
O(1,3) = (K1^2*K2^2+X16o^2*(K2+X4o)^2+K1*K2*(K2*(WNT0+2*X16o)+2*(WNT0+X16o)*X4o))/(K1*K2+X16o*(K2+X4o))^2 ;

O(2,1) = 1 ;
O(2,2) = 1 ;
O(2,3) = 0 ;

O(3,1) = -1 ;
O(3,2) = (K2*WNT0*X16o*(K1+X16o))/(K1*K2+X16o*(K2+X4o))^2 ;
O(3,3) = (K1*K2*WNT0*X4o)/(K1*K2+X16o*(K2+X4o))^2 ;

d = zeros(3,1);
d(1,1) = k17*X14o*X15o - k18*X16o - k4*X4o + k5*X5o + k3*X6o*X7o - km3*X4o ;
d(2,1) = -dCdt0(X5) - dCdt0(X6) - dCdt0(X9) ;
d(3,1) = 0 ;

z = O\d;

dCdt0(X3) = z(1);
dCdt0(X4) = z(2);
dCdt0(X16) = z(3);

% the rest (dependent variables):

dCdt0(X1) = (-K1*K2*WNT0*X16o/(K1*K2+X16o*(K2+X4o))^2)*dCdt0(X4) + (-K1*K2*WNT0*(K2+X4o)/(K1*K2+X16o*(K2+X4o))^2)*dCdt0(X16) ;
dCdt0(X2) = (-K2*WNT0*X16o^2/(K1*K2+X16o*(K2+X4o))^2)*dCdt0(X4) + (K1*K2^2*WNT0/(K1*K2+X16o*(K2+X4o))^2)*dCdt0(X16) ;

%fixed_dCdt0 = [];
fixed_dCdt0 = zeros(1,var);

fixed_dCdt0(X1) = 0;
fixed_dCdt0(X2) = 0;
fixed_dCdt0(X3) = 0;
fixed_dCdt0(X4) = 0;
fixed_dCdt0(X5) = 0;
fixed_dCdt0(X6) = 0;
fixed_dCdt0(X7) = 0;
fixed_dCdt0(X8) = 0;
fixed_dCdt0(X9) = 0;
fixed_dCdt0(X10) = 0;
fixed_dCdt0(X11) = 0;
fixed_dCdt0(X12) = 0;
fixed_dCdt0(X13) = 0;
fixed_dCdt0(X14) = 0;
fixed_dCdt0(X15) = 0;
fixed_dCdt0(X16) = 0;
fixed_dCdt0(X17) = 0;
fixed_dCdt0(X18) = 0;
fixed_dCdt0(X19) = 0;
fixed_dCdt0(SigER) = 0;
fixed_dCdt0(SigM) = 0;
fixed_dCdt0(SigERC) = 0;

% To get consistent C and dCdt initial conditions
[C0mod,dCdt0mod] = decic(@rcn_functions,tspan,C0,fixed_C0,dCdt0,fixed_dCdt0);

% Solve system of DAEs
% [t,S] = ode15i(@rcn_functions,tspan,C0mod,dCdt0mod);

% For cases where ode15i fails (change error tolerance):
options=odeset('AbsTol',1e-6,'RelTol',1e-3);

[t,S] = ode15i(@rcn_functions,tspan,C0mod,dCdt0mod,options);

conc_out = S(tr,:);



function daes = rcn_functions(t,C,dCdt)

% for scalar t and column vector C and dCdt
% returns column vector corresponding to daes(t,C,CP)
% (daes = differential algebraic equations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RCN Parameters
% found in m-files: initialize_rcn.m, solve_rcn.m, rcn_functions
% (withing solve_rcn)

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

X1 = 1;     X2 = 2;     X3 = 3;     X4 = 4;     X5 = 5;     
X6 = 6;     X7 = 7;     X8 = 8;     X9 = 9;     X10 = 10;
X11 = 11;   X12 = 12;   X13 = 13;   X14 = 14;   X15 = 15;   
X16 = 16;   X17 = 17;   X18 = 18;   X19 = 19;
SigER = 20;     SigM = 21;      SigERC = 22;


daes = zeros(size(C));


daes(1) = WNT0 - C(X1) - C(X2) - C(X3) ;

daes(2) = WNT0*K2*C(X16)/(K1*K2 + K2*C(X16) + C(X4)*C(X16)) - C(X2) ;

daes(3) = WNT0*C(X4)*C(X16)/(K1*K2 + K2*C(X16) + C(X4)*C(X16)) - C(X3) ;

daes(4) = km3*C(X4) - k3*C(X6)*C(X7) - dCdt(X6) ;

daes(5) = APC0/(1+C(X10)/K8) - C(X7) ;

daes(6) = C(X7)*C(X10)/K8 - C(X8) ;

daes(7) = C(X5)*C(X10)/K6 - C(X9) ;

daes(8) = K11*TCF0/(K11+C(X10)) - C(X11) ;

daes(9) = TCF0*C(X10)/(K11+C(X10)) - C(X12) ;

daes(10) = Tmax*(C(X12)^yps/(Krna^yps + C(X12)^yps)) - k13*C(X13) - dCdt(X13) ;

daes(11) = Pmax*C(X13) - k19*C(X14) - dCdt(X14) ;

daes(12) = v15 - k17*C(X14)*C(X15) - k16*C(X15) - dCdt(X15) ;

daes(13) = v20 - k21*C(X17) - dCdt(X17) ;

daes(14) = k21*C(X17) - k22*C(X18) + k23*C(X19) - dCdt(X18) ;

daes(15) = k22*C(X18) - k23*C(X19) - k25*C(X19) - dCdt(X19) ;

daes(16) = (v20/C(X17))*((1 - VtG*C(X14)/(Km+C(X14))) - C(SigER)) - dCdt(SigER) ;

daes(17) = (k21*C(X17)/C(X18))*(C(SigER)-C(SigM)) + (k23*C(X19)/C(X18))*(C(SigERC)-C(SigM)) - dCdt(SigM) ;

daes(18) = (k22*C(X18)/C(X19))*(C(SigM)-C(SigERC)) - dCdt(SigERC) ;

daes(19) = -k4*C(X4) + k5*C(X5) - k7*C(X9) + v9 + k25*C(X19) - k10 + dCdt(X5) - (APC0*K8/(K8+C(X10))^2 + K11*TCF0/(K11+C(X10))^2 + 1)*dCdt(X10) ;

daes(20) = k4*C(X4) - k5*C(X5) - (C(X10)/K6 + 1)*dCdt(X5) - (C(X5)/K6)*dCdt(X10) ;

daes(21) = k17*C(X14)*C(X15) - k18*C(X16) - k4*C(X4) + k5*C(X5) + k3*C(X6)*C(X7) - km3*C(X4)...
    - (K1*K2*WNT0*C(X16)/(K1*K2+C(X16)*(K2+C(X4)))^2 + K2*WNT0*C(X16)*(K1+C(X16))/(K1*K2+C(X16)*(K2+C(X4)))^2 + 1)*dCdt(X4)...
    - (K1*K2*WNT0*(K2+C(X4))/(K1*K2+C(X16)*(K2+C(X4)))^2 + K1*K2*WNT0*C(X4)/(K1*K2+C(X16)*(K2+C(X4)))^2 + 1)*dCdt(X16) ;

daes(22) = dCdt(X3) + dCdt(X4) + dCdt(X5) + dCdt(X6) + dCdt(X9);

