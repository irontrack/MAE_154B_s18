% CONVERT VN DIAGRAM DATA TO ALPHA 
% THE FOLLOWING PROGRAM OUTPUTS AN ARRAY CONTAINING
% APLHA, CL, CD, AND CM FOR ALL CRITICAL POINTS
% THE FIRST 6 ROWS CONTAIN POINT DATA FOR SEA LEVEL
% THE LAST 6 CONTAIN POINT DATA FOR 12K FEET.

% HERE IS A TEMPLATE

%   NOTE ****ALPHA*** IS IN ****RADIANS ***

%           ALPHA        CL       CD         CM
% PHAA SEA
% PGAA SEA
% PLAA SEA
% NLAA SEA
% NGAA SEA
% NHAA SEA
% REPEAT IN SAME ORDER FOR 12K FEET
% A DIAGRAM OF ACCOCIATED POINTS IS IN A .JPG NAMMED "CRITICAL POINTS WITH LABEL" 
function alpha = vn_find_alpha()
    
    % note these velocities are in mps and need to be converted to fps
    % later.
    vn_points = [50.4361, 4.4;...
                87.5, 6.9559;...
                131.2, 5.4849;...
                131.2, -4.4849;...
                87.5, -5.9559;...
                40.3401, -1.76;...
                61.6608, 4.4;...
                87.5, 6.9477;...% to do
                131.2, 5.4795;... % to do
                131.2, -4.4795;...% to do
                87.5, -5.9477;...% to do
                49.3673, -1.76];
                
    b = 37;
    w = 2000/b;     %[pounds/ft]
    c = 4.41667;  %[ft]
    e = 0.79;
    AR = b/c;
    
    CL_alpha_12k = 4.8716;
    CL_0_12k = 0.2364;
    CL_alpha_sea = 5.0646;
    CL_0_sea = 0.2445;
    CD_0 = 0.02;
    
    % conversion factor feet to meeters
    ftm = 3.0480e-01;
    
    % conversion factor rho_si to rho_imp
    
    rho_imp = (2.20462/(32.174*35.3147));
    
    rho_sea = rho_imp*density(0);
    rho_12k = rho_imp*density(12000*ftm);
    
    temp = zeros(12,4);
    % convert meters per second to feet per second
    mps_fps = 3.28084;
    for i = 1:12
        
        v = mps_fps*vn_points(i,1);
        n = vn_points(i,2);
        
        L = w*n;
        if i <7
            q = 0.5*rho_sea*v^2;
            temp(i,1) = ( L/(c*q) - CL_0_sea)/CL_alpha_sea;
            temp(i,2) = CL_alpha_sea*temp(i,1) + CL_0_sea;
            temp(i,3) = CD_0 + (1/(pi*AR*e))*(temp(i,2)^2);
            
        else
            q = 0.5*rho_12k*v^2;
            temp(i,1) = ( L/(c*q) - CL_0_12k)/CL_alpha_12k;
            temp(i,2) = CL_alpha_12k*temp(i,1) + CL_0_12k;
            temp(i,3) = CD_0 + (1/(pi*AR*e))*(temp(i,2)^2);
        end
        
    end
        for i = 1:12
            temp(i,4) = 0.056;
        end
        alpha = temp;
        
end
    
    
    