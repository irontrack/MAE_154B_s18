% this script prints out data for 3d lift curve slope of the 
% naca2412 airfoil at sealevel and 120000 ft.


% table of data on airfoil at 12000ft.
maxalt = naca_table_maxalt();

% table of data on airfoil at sealvl.
sealvl = naca_table_sealvl();

% vector of alpha values from Xfoil data

% note: only rows 18-54 are taken because this is the gaurenteed linear
% regime.

alpha_max = maxalt(18:54,1);
% vector of corresponding Cl values from Xfoil data
Cl_max = maxalt(18:54,2);

% vector of alpha values from Xfoil data
alpha_sea = sealvl(18:54,1);
% vector of corresponding Cl values from Xfoil data
Cl_sea = sealvl(18:54,2);

%getting 2d liftcurve slope from data
temp = polyfit(alpha_max,Cl_max,1);
Cl_alpha_max = (180/pi)*temp(1);
Cl_0_max = temp(2);

temp = polyfit(alpha_sea,Cl_sea,1);
Cl_alpha_sea = (180/pi)*temp(1);
Cl_0_sea = temp(2);


%calculate 3d liftcurve slope data

A = 8.3774;
e = 0.79;

CL_alpha_max = Cl_alpha_max/(1+(Cl_alpha_max/(pi*e*A)));
CL_alpha_sea = Cl_alpha_sea/(1+(Cl_alpha_sea/(pi*e*A)));

alpha_max = (pi/180)*maxalt(:,1);
alpha_sea = (pi/180)*sealvl(:,1);
CL_max = CL_alpha_max.*alpha_max + Cl_0_max;
CL_sea = CL_alpha_sea.*alpha_sea + Cl_0_sea;

CL_pmax = CL_alpha_max*(16)*(pi/180) + Cl_0_max;
CL_nmax = CL_alpha_max*(-14.5)*(pi/180) + Cl_0_max;

% plot CL vs alpha at sea level
plot(alpha_max,CL_max,'r')
hold on
plot(alpha_max,maxalt(:,2))


%display final results
disp('altitude    CL_alpha    CL_0    ');
fprintf('\n12000ft     %8.4f',CL_alpha_max);
fprintf('    %8.4f',Cl_0_max);
fprintf('\nsea lvl     %8.4f',CL_alpha_sea);
fprintf('    %8.4f\n',Cl_0_sea);
fprintf('CL_pmax %8.4f\n',CL_pmax);
fprintf('CL_nmax %8.4f\n',CL_nmax);


