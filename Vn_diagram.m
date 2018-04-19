%Script for constructing V-n diagram.  
clc;
close all;

%define aero constants (you should check for the correctness)
W = 2000;
S = 163;
rho = 0.00238; %@ sea level or 12k ft
Vc = 170;
Vd = 255;
AR = 8.3774;
e=0.79;
CD0 = 0.00549;
c=4.4167;
g=32.2;
b = 18.5;  %semi-span

%%% ToDo Task:
%    %% Plots for 2D C_l and C_d vs. alpha... 
     %...Use your own airfoil data from XFoil for plotting to determine 2D
     % lift and drag curve slopes.
%Positive alphas
% figure 
% plot(alpha1, Cl1);
% title('Cl1 vs. Alpha');
% 
% figure
% plot(alpha1, Cd1);
% title('Cd1 vs. Alpha');
% 
% %Negative alphas
% figure 
% plot(alpha2, Cl2);
% title('Cl2 vs. Alpha');
% 
% figure
% plot(alpha2, Cd2);
% title('Cd2 vs. Alpha');
% %%%%%%%%%%%%%
%%%


%After you've determined the 2D lift and drag slopes from the above plots,...
%...convert them to 3D slopes. The following numerical values used are 
%arbitrarily chosen, so use your own numbers.

%Positive alphas
amax_p = 17.5/180*pi;
% Cl_alpha_p = 0.103*180/pi; %2D slope
% CL_alpha_p = Cl_alpha_p/(1+Cl_alpha_p/(pi*AR*e)); %3D
CL_alpha_p = 4.8719;
CL_p_max = CL_alpha_p*amax_p;
CD_p_max = CD0 + CL_p_max^2/(pi*AR*e);
Cz_max = CL_p_max*cos(amax_p)+CD_p_max*sin(amax_p); %3D

mu = 2*W/(rho*g*c*CL_alpha_p*S);
kg = 0.88*mu/(5.3+mu);   %gust alleviation factor

%Negative alphas
amax_n = -17/180*pi;
% Cl_alpha_n = 0.11289*180/pi;
% CL_alpha_n = Cl_alpha_n/(1+Cl_alpha_n/(pi*AR*e));
CL_alpha_n = 4.8719;
CL_n_max = CL_alpha_n*amax_n;
CD_n_max = CD0 + CL_n_max^2/(pi*AR*e);
Cz_neg = CL_n_max*cos(amax_n)+CD_n_max*sin(amax_n);

%Maneuver Limit	Envelope load factors
v = 0:270; %mph
vfps = v*5280/3600; %fps
n_lift = 1.25*(0.5*rho*S*Cz_max/W).*vfps.^2;
n_neg =  (0.5*rho*S*Cz_neg/W).*vfps.^2;


lim_load_plus = 4.4;
lim_load_minus = -1.76;


%Gust Limit Envelope
g_50_plus = 1+kg*rho*50*CL_alpha_p*S/(2*W).*vfps;
g_25_plus = 1+kg*rho*25*CL_alpha_p*S/(2*W).*vfps;
g_50_minus = 1-kg*rho*50*CL_alpha_p*S/(2*W).*vfps;
g_25_minus = 1-kg*rho*25*CL_alpha_p*S/(2*W).*vfps;


%%% ToDo Task:
%Define the positive and negative stall speeds
%%%

%%% ToDo Task:
%Limit	Combined  Envelope... 
%Include (i.e. graph) the Allowable Envelope in the VN diagrams.
%Total allowable envelope for sea level/12k ft is obtained by adding 
%the maneuvering load and gust load limits. 		
%%%


figure; grid on; hold on;set(gcf,'color',[1 1 1]);


tempInd = find(n_lift<lim_load_plus); plot(v(tempInd),n_lift(tempInd),'linewidth',2); 
tempInd = find(n_neg>lim_load_minus); plot(v(tempInd),n_neg(tempInd),'linewidth',2);

tempInd = find(n_lift>lim_load_plus,1); plot([v(tempInd) v(end)],[lim_load_plus lim_load_plus],'linewidth',2);

tempInd = find(n_neg<lim_load_minus,1); plot([v(tempInd) Vc],[lim_load_minus lim_load_minus],'linewidth',2);

plot([Vd Vd],[-1 lim_load_plus],'linewidth',2)
plot([Vc Vd],[lim_load_minus -1],'linewidth',2)
plot(v(1:230),g_50_plus(1:230),'--g','linewidth',2)
plot(v(1:270),g_25_plus(1:270),'--g','linewidth',2)
plot(v(1:230),g_50_minus(1:230),'--g','linewidth',2)
plot(v(1:270),g_25_minus(1:270),'--g','linewidth',2)
plot([Vc Vd],[g_50_plus(230) g_25_plus(270)],'--g','linewidth',2)
plot([Vc Vd],[g_50_minus(230) g_25_minus(270)],'--g','linewidth',2)




plot([230 230],[lim_load_minus lim_load_plus],'--r','linewidth',1)
xlabel('V (mph)','fontsize',16,'fontweight','bold');ylabel('n','fontsize',16,'fontweight','bold')
set(gca,'FontSize',16,'fontweight','bold');
ylim([-3 6])
title('V-N Diagram for NACA 2412');
print('-dpng',[pwd,'figure',num2str(1),'.png']);

%%% ToDo Task:
%construct a VN diagram for both sea level and 12k ft
%%%


