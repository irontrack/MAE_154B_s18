CL_pmax_12k =  1.5969;
CL_nmax_12k = -0.9965;
CL_pmax_sea = 1.6588;
CL_nmax_sea = -1.0372;

 b = 37*3.0480e-01;       %[m]
 w = (2000/2.2)*9.81;   %[N]
 c = 4.41667*3.0480e-01;  %[m]
 s = b*c;
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
            
 V_1sea = linspace(0,vn_points(1,1),100);
 V_112k = linspace(0,vn_points(7,1),100);
 V_6sea = linspace(0,vn_points(6,1),100);
 V_612k = linspace(0,vn_points(12,1),100);
 
 
 rho_sea = 1.2251;  % [kg/m^3]
 rho_12k = 0.8514;  % [kg/m^3]
 n1_sea = (0.5*rho_sea.*V_1sea.^2)*s*CL_pmax_sea/(w);
 n2_sea = [4.4 4.4 -1 -1.76 -1.76];
 n6_sea = (0.5*rho_sea.*V_6sea.^2)*s*CL_nmax_sea/w;
 V_2sea = [vn_points(1,1) 131.2 131.2 87.5 vn_points(6,1)];
 
 n1_12k = (0.5*rho_12k.*V_112k.^2)*s*CL_pmax_12k/w;
 n2_12k = [4.4 4.4 -1 -1.76 -1.76];
 n6_12k = (0.5*rho_12k.*V_612k.^2)*CL_nmax_12k*s/w;
 V_212k = [vn_points(7,1) 131.2 131.2 87.5 vn_points(12,1)];
 
 n_gust_sea = [1 vn_points(2,2) vn_points(3,2) 1 vn_points(4,2) vn_points(5,2) 1];
 v_gust_sea = [1 vn_points(2,1) vn_points(3,1) 1 vn_points(4,1) vn_points(5,1) 1];
 n_gust_12k = [1 vn_points(8,2) vn_points(9,2) 1 vn_points(10,2) vn_points(11,2) 1];
 v_gust_12k = [1 vn_points(8,1) vn_points(9,1) 1 vn_points(10,1) vn_points(11,1) 1];
 
 plot(V_1sea,n1_sea,'b',V_2sea,n2_sea,'b',V_6sea,n6_sea,'b')
 hold on
 plot(v_gust_sea,n_gust_sea,'g--')

 grid on
xlabel('velocity (m/s)')
ylabel('load')
title('at sea level')
 figure (2)

 plot(V_112k,n1_12k,'b',V_212k,n2_12k,'b',V_612k,n6_12k,'b')
  hold on
 plot(v_gust_12k,n_gust_12k,'g--')
 
 grid on
 xlabel('velocity (m/s)')
ylabel('load')
title('at 12k feet')

            