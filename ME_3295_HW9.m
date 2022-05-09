%% ME 3295 HOMEWORK 9
%
% Student:          Kevin Cayo
% Instructer:       Alexei Poludnenko
% Courese:          ME 3295
% Date:             Friday, April 23, 2022
%%

%% Shock Mach Number

clc
clear

%% Knowns
gamma = 1.4;
R = 287; % j/(kg*K)

%% Station 3 Pre shock
P3 = 1   % atm
T3 = 273 % k
v3 = 0
a3 = sqrt(gamma*R*T3)

%% Station 2 Post Shock
P2 = 4.5 %atm

%% Station 1 Throat
M1 = 1


%% station 0 Reservoir
v0 = 0

%% Pressure ratios
P2_P3 = P2/P3 % 4.5, go to appendix B, Normal shock Relations
M3 = 2          % pre
M2 = .5774      % post
T2_T3 = 1.687   % Post/Pre
T2 = T3 * T2_T3

%% finding Stagnation Pressure & Temperature at station 2 & 3
T02 = T2*(1+((gamma-1) / 2) *M2^2)
P02 = P2*(1+((gamma-1) / 2) *M2^2)^((gamma)/(gamma-1))
T03 = T3*(1+((gamma-1) / 2) *M3^2)
P03 = P3*(1+((gamma-1) / 2) *M3^2)^((gamma)/(gamma-1))

% print(T02,T03,P02, P03)

%% Sonic Conditons at station 3-2
% - For gamma = 1.4 -
Tstar_T0 = (1+((gamma-1) / 2) *1^2)^(-1)
Pstar_P0 = 1/(1+((gamma-1) / 2) *1^2)^((gamma)/(gamma-1))
Tstar3 = Tstar_T0 * T03
Pstar3 = Pstar_P0 * P03

Tstar2 = Tstar_T0 * T02
Pstar2 = Pstar_P0 * P02

Mstar3 = sqrt(((gamma+1)*M3^2)/(2+(gamma-1)*M3^2))
Mstar2 = sqrt(((gamma+1)*M2^2)/(2+(gamma-1)*M2^2))

%% Calculating Velocities
Tstar1 = Tstar2
a1 = sqrt(gamma*R*Tstar1);
a2 = sqrt(gamma*R*T2);
u3 = M3 * a3;
u2 = M2 * a2;
u1 = M1 * a1;
V2 = u3 - u2
V1 = u3 - u1





