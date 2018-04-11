function dat = naca_table_maxalt()
    % % creates a table with max altitude (12000ft) data on the naca 
    % 2412 airfoil.
    
    % @ param table quatities include Cl, Cd, alpha
    % Cdp Cm top_Xtr, and bot_xtr
            %alpha     Cl          Cd          Cdp        Cm
            
            %                                                           
    naca2412_maxalt = fopen('naca2412_maxalt.pol');
    A = fscanf(naca2412_maxalt,'%f',[7 Inf]);
    fclose(naca2412_maxalt);
    
    dat = transpose(A);
end
