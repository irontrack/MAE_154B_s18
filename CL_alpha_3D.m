% this script prints out data for 3d lift curve slope of the 
% naca2412 airfoil at sealevel and 120000 ft.


% table of data on airfoil at 12000ft.
maxalt = naca_table_maxalt();

% table of data on airfoil at sealvl.
sealvl = naca_table_sealvl();

% vector of alpha values from Xfoil data
alpha_max = maxalt(:,1);
% vector of corresponding Cl values from Xfoil data
Cl_max = maxalt(:,2);

% vector of alpha values from Xfoil data
alpha_sea = sealvl(:,1);
% vector of corresponding Cl values from Xfoil data
Cl_sea = sealvl(:,2);

%getting 2d liftcurve slope from data
temp = polyfit(alpha_max,Cl_max,1);
Cl_alpha_max = temp(1);
Cl_0_max = temp(2);

temp = polyfit(alpha_sea,Cl_sea,1);
Cl_alpha_sea = temp(1);
Cl_0_sea = temp(2);


%calculate 3d liftcurve slope data

A = 8.3774;
e = 0.79;

CL_alpha_max = Cl_alpha_max/(1+(Cl_alpha_max/(pi*e*A)));
CL_alpha_sea = Cl_alpha_sea/(1+(Cl_alpha_sea/(pi*e*A)));

disp('altitude    CL_alpha    CL_0    ');
fprintf('\n12000ft     %8.4f',CL_alpha_max);
fprintf('    %8.4f',Cl_0_max);
fprintf('\nsea lvl     %8.4f',CL_alpha_sea);
fprintf('    %8.4f\n',Cl_0_sea);

