%% ME 3295 HOMEWORK 6
%
% Student:          Kevin Cayo
% Instructer:       Alexei Poludnenko
% Courese:          ME 3295
% Date:             Sunday, March 20, 2022
%%
clc
clear

%% PARAMETERS / KNOWN VARIABLES / EQUATIONS

rho = 0.3643;           %kg/m^3
v = 200 * 0.44704;      %Velocity, m/s
c = 1.2;                %Chord, m
h = 5000;               %altitude,

% Airfoil shape
m = 0.02;
p = 0.4;
t = 0.12;

%% EQUATIONS
zc_f = @(x) m/p.^2 * (2*p*x - x.^2) ;
zc_r = @(x) m/ (1-p) ^2 *((1-2*p) + 2*p*x -x.^2) ;
dzc_dx_f = @(x) 2 * m / p^2 *(p-x);
dzc_dx_r = @(x) 2 * m / (1-p)^2 * (p-x);
x = @(theta) 0.5 * (1- cos(theta));
theta = linspace(0, pi() , 1000) ;
theta_p = acos(1-2*p) ;
i_mid = find(theta>=theta_p,1, 'first');

%% ZERO LIFT ANGLE OF SHOCK
irng = 1:i_mid-1;
integrand(irng) = dzc_dx_f(x(theta(irng))).* (cos(theta(irng))-1) ;
irng = i_mid: length(theta) ;
integrand(irng) = dzc_dx_r (x(theta (irng))).* (cos(theta(irng))-1) ;
a10 = -1/pi ()*trapz(theta, integrand) ;
a10_d = rad2deg(a10) ;

cl = @(a) 2*pi () * (a-al0) ;

% Plot cl
arng = linspace (-10,10, 1000) ;
plot(arng,cl(deg2rad(arng)), 'b') ;
grid on
xlabel ('\alpha (deg)')
ylabel ('c_l')
title ('Lift Coefficient - NACA 2412')

%% CENTER PRESURE
% A1
irng = 1:i_mid-1;
integrand (irng) = dzc_dx_f(x(theta(irng))).*(cos(theta(irng))) ;
irng = i_mid:lenath(theta)
integrand (irng) = dzc_dx_r(x(theta(irng))).*(cos(theta(irng))) ;
A1 = 2/pi()*trapz(theta,integrand);

%A2
irng = 1:i_mid-1;
integrand(irng) = dzc_dx_f(x(theta(irng))).*(cos(2*theta(irng))) ;
irng = i_mid:length (theta) ;
integrand (irng) = dzc_dx_r(x(theta(irng))).*(cos(2*theta(irng))) ;
A2= 2/pi()*trapz (theta,integrand);

x_cp = c/4 * (1 + pi()./cl(deg2rad(arng))*(A1-A2)) ;
figure
plot(arng, x_cp, 'b')
hold on
plot (arng,ones(size(arng)) *c/4, 'b--')
grid on
xlabel (' \alpha (deg) ')
label ('x_{cp} (m) ')
title ('Center of Pressure - NACA 2412')
ylim([-1 1])

m_c4 = pi()/4*(A2-A1) ;