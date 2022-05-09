%% ME 3295 HOMEWORK 10
%
% Student:          Kevin Cayo
% Instructer:       Alexei Poludnenko
% Courese:          ME 3295
% Date:             Friday, April 29, 2022
%%
clc
clear

%% PARAMETERS / KNOWN VARIABLES / EQUATIONS

M1=8.0; %Change to 1.5 for part b
gam=1.4;
[max_alpha_oblique, ~] = max_alpha_oblique_shock(M1,gam);
max_alpha_expanfan = max_alpha_expan_fan(M1,gam);
limiting_alpha = min(max_alpha_oblique,max_alpha_expanfan);
alpha_deg = 0:0.1:limiting_alpha;
alpha = deg2rad(alpha_deg);

%Expansion Fan
[~,nu_M1,~] = flowprandtlmeyer(gam,M1);
nu_M2 = alpha_deg + nu_M1;
M2 = zeros(size(nu_M2));
for ind = 1:length(nu_M2)
 M2(ind) = flowprandtlmeyer(gam,nu_M2(ind),'nu');
end
Pt2_P2 = (1+((gam-1)/2)*(M2.^2)).^(gam/(gam-1));
Pt1_P1 = (1+((gam-1)/2)*(M1^2)).^(gam/(gam-1));
P2_P1 = (1./Pt2_P2).*Pt1_P1;

%Oblique Shock
beta = oblique_beta(M1,alpha_deg,gam,0);
Mn1 = M1*sin(deg2rad(beta));
P3_P1 = 1+(((2*gam)/(gam+1)).*((Mn1.^2)-1));

%% GENERATING PLOTS
% Static Pressure Ratios vs Angle of Attack
figure(1)
hold on
plot(alpha_deg, P3_P1)
ylabel('Static Pressure ratio')
yyaxis right
plot(alpha_deg,P2_P1)
hold off
title('Static Pressure Ratios vs Angle of Attack')
legend('P3/P1','P2/P1')
xlabel('Angle of Attack [degrees]')
ylabel('Static Pressure Ratio')

% Aerodynamic Coefficients
figure(2)
Cl = (2/(gam*(M1^2))).*(P3_P1-P2_P1).*cos(alpha);
Cd = (2/(gam*(M1^2))).*(P3_P1-P2_P1).*sin(alpha);
figure(2)
hold on
plot(alpha_deg,Cd)
plot(alpha_deg,Cl)
ylabel('Coefficient')
yyaxis right
plot(alpha_deg,Cl./Cd)
ylabel('C/D Ratio')
ylim([0 25])
hold off
title('Aerodynamic Coefficients')
legend('Cd','Cl','L/D')
xlabel('Angle of Attack (degrees)')

function Beta=oblique_beta(M,theta,gamma,n)
 Beta=zeros(size(theta));
 for ind = 1:length(theta)
 theta(ind)=theta(ind)*pi/180; 
 mu=asin(1/M); 
 c=tan(mu)^2;
 a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta(ind));
 b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta(ind));
 d=sqrt(4*(1-3*a*b)^3/((27*a^2*c+9*a*b-2)^2)-1);
 Beta(1,ind)=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a^2*c+9*a*b-2))...
 /(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)))*180/pi;
 end
end

function [max_alpha,beta_at_alpha] = max_alpha_oblique_shock(M1,gamma)
 gam=gamma;

 beta_deg = 0.1:0.1:90;
 beta = deg2rad(beta_deg);

 theta = acot(tan(beta).*((((gam+1)*(M1^2))./(2*((M1^2)*((sin(beta).^2))-1)))-1));
 theta_deg=rad2deg(theta);

 [max_alpha,index] = max(theta_deg);
 beta_at_alpha = beta_deg(index);
end

function max_alpha = max_alpha_expan_fan(M1,gamma)
 gam=gamma;
 v_M = rad2deg((((gam + 1)/(gam - 1))^(1/2))...
 *(atan((((gam - 1)/(gam + 1))*((M1^2) - 1))^(1/2)))...
 - atan(((M1^2) - 1)^(1/2)));
 v_max = rad2deg((pi/2)*((sqrt((gam+1)/(gam-1)))-1));
 max_alpha = v_max-v_M;
end
