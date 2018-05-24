% gust load factor calculator
w = 2000;
c = 4.41667;    %[ft]
b = 37; %[ft]
s = b*c;    %ft^2]
g = 32.174;  %[ft/s^2]
Vc = 87.5;   %[m/s]
Vd = 131.02; %[m/s]

%conversion factor from [m/s] to knotts
m2k = 1.94384; 

Uec = 50;   %[ft/s]
Ued = 25;   %[ft/s]

alpha_sea = 5.0646; %[/rad]
alpha_12k = 4.8716; %[/rad]


% conversion factor rho_si to rho_imp
rho_imp = (2.20462/(32.174*35.3147));
rho_sea = 1.2251*rho_imp;  % [slug/ft^3]
rho_12k = 0.8514*rho_imp;  % [slug/ft^3]

mu_sea = 2*(w/s)/(rho_sea*c*alpha_sea*g);
mu_12k = 2*(w/s)/(rho_12k*c*alpha_12k*g);

kg_sea = 0.88*mu_sea/(5.3 +mu_sea);
kg_12k = 0.88*mu_12k/(5.3 +mu_12k);

% positive gust load factors
ngc_sea = 1 + kg_sea*alpha_sea*Uec*Vc*m2k/(498*(w/s));
ngd_sea = 1 + kg_sea*alpha_sea*Ued*Vd*m2k/(498*(w/s));

ngc_12k = 1 + kg_sea*alpha_12k*Uec*Vc*m2k/(498*(w/s));
ngd_12k = 1 + kg_sea*alpha_12k*Ued*Vd*m2k/(498*(w/s));

vn = [0 1;
      87.5 ngc_sea;
      131.2 ngd_sea;
      131.2 (1 - ngd_sea);
      87.5 (1 - ngc_sea);
      0 1;
      0 1;
      87.5 ngc_12k;
      131.2 ngd_12k;
      131.2 (1 - ngd_12k);
      87.5 (1 - ngc_12k);
      0 1];
      