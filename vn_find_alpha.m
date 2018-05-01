% CONVERT VN DIAGRAM DATA TO ALPHA 

function alpha = vn_find_alpha()

    vn_data = fopen('vn_data.txt');
    vn_points = fscanf(vn_data,'%f',[2,inf]);
    fclose(vn_data);
    b = 37;
    w = 2000/b;     %[pounds/ft]
    c = 4.41667;  %[ft]
    
    
    CL_alpha_12k = 4.8716;
    CL_0_12k = 0.2364;
    CL_alpha_sea = 5.0646;
    CL_0_sea = 0.2445;
    
    % conversion factor feet to meeters
    ftm = 3.0480e-01;
    
    % conversion factor rho_si to rho_imp
    
    rho_imp = (2.20462/(32.174*35.3147));
    
    rho_sea = rho_imp*density(0);
    rho_12k = rho_imp*density(12000*ftm);
    
    temp = zeros(12,1);
    
    for i = 1:12
        
        v = vn_points(1,i);
        n = vn_points(2,i);
        
        L = w*n;
        if i <7
            q = 0.5*rho_sea*v^2;
            temp(i) = ( L/(c*q) - CL_0_sea)/CL_alpha_sea;
        else
            q = 0.5*rho_12k*v^2;
            temp(i) = ( L/(c*q) - CL_0_sea)/CL_alpha_12k;
        end
    end
        alpha = temp;
        
end
    
    
    