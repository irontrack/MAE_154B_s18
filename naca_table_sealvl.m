function dat = naca_table_sealvl()
    % creates a table with sealvl data on the naca 
    % 2412 airfoil.
    
    % @ param table quatities include Cl, Cd, alpha
    % Cdp Cm top_Xtr, and bot_xtr
    
    
    naca2412_sealvl = fopen('naca2412_sealvl.pol');
    A = fscanf(naca2412_sealvl,'%f',[7 Inf]);
    fclose(naca2412_sealvl);
    
    dat = transpose(A);
    
    
end