CL_pmax_12k =  1.5969;
CL_nmax_12k = -0.9965;
CL_pmax_sea = 1.6588;
CL_nmax_sea = -1.0372;

 rho_sea = 1.2251;  % [kg/m^3]
 rho_12k = 0.8514;  % [kg/m^3]
 b = 37*3.0480e-01;       %[m]
 w = (2000/2.2)*9.81;   %[N]
 c = 4.41667*3.0480e-01;  %[m]
 s = b*c;
 
Vsp_sea = sqrt(2*w/(rho_sea*CL_pmax_sea*s));
Vsn_sea = sqrt(-2*w/(rho_sea*CL_nmax_sea*s));
Vsp_12k = sqrt(2*w/(rho_12k*CL_pmax_12k*s));
Vsn_12k = sqrt(-2*w/(rho_12k*CL_nmax_12k*s));

% no gust points
vn_points =[50.4361, 4.4;...
            131.2, 4.4;...
            131.2, -1;...
            87.5, -1.76;...
            40.3401, -1.76;...
            61.6608, 4.4;...
            131.2, 4.4;...
            131.2, -1;...
            87.5, -1.76;...
            49.3673, -1.76];
        
% flight envelope points
vn_envelope = [65.5623, 4.4;...
                87.5000, 5.5378;...
                131.2, 4.4;...
                131.2000,-3.3974;...
                87.5000, -4.5378;...
                43.6094, -1.76;...
                68.1568, 4.4;...
                87.5000,5.3649;...% to do
                131.2, 4.4;... % to do
                131.2000, -3.2679;...% to do
                87.5000, -4.3649;...% to do
                49.3673, -2.0269];
 vn_envelope = transpose(vn_envelope);
% gust points            
vn_gust =      [ 0    1.0000;...
               87.5000    5.5378;...
              131.2000    4.3974;...
              131.2000   -3.3974;...
               87.5000   -4.5378;...
                     0    1.0000;...
                     0    1.0000;...
               87.5000    5.3649;...
              131.2000    4.2679;...
              131.2000   -3.2679;...
               87.5000   -4.3649;...
                     0    1.0000]; 
                 
                 
 V_1sea = linspace(0,vn_points(1,1),100);
 V_112k = linspace(0,vn_points(6,1),100);
 V_6sea = linspace(vn_points(5,1),0,100);
 V_612k = linspace(vn_points(10,1),0,100);
 % find indicies of stall angle
 i1sea = find(abs(V_1sea - Vsp_sea)<.4);
 i6sea = find(abs(V_6sea - Vsn_sea)<.4);
 i112k = find(abs(V_112k - Vsp_12k)<.4);
 i612k = find(abs(V_612k - Vsn_12k)<.4);
 
 n1_sea = (0.5*rho_sea.*V_1sea.^2)*s*CL_pmax_sea/(w);
 n6_sea = (0.5*rho_sea.*V_6sea.^2)*s*CL_nmax_sea/w;
 n_points_sea = [n1_sea, vn_points(2,2), vn_points(3,2), vn_points(4,2), vn_points(5,2), n6_sea];
 v_points_sea = [V_1sea, vn_points(2,1), vn_points(3,1), vn_points(4,1), vn_points(5,1), V_6sea];

 
 
 n1_12k = (0.5*rho_12k.*V_112k.^2)*s*CL_pmax_12k/w;
 n6_12k = (0.5*rho_12k.*V_612k.^2)*s*CL_nmax_12k/w;
 n_points_12k = [n1_12k, vn_points(7,2), vn_points(8,2), vn_points(9,2), vn_points(10,2), n6_12k];
 v_points_12k = [V_112k, vn_points(7,1), vn_points(8,1), vn_points(9,1), vn_points(10,1), V_612k];
 
 % gust plots
 
 ng1_sea = [0 n1_sea(i1sea:length(n1_sea)) vn_envelope(2,1:6) n6_sea(1:i6sea) 0 0];
 vg1_sea = [Vsp_sea V_1sea(i1sea:length(V_1sea)) vn_envelope(1,1:6) V_6sea(1:i6sea) Vsn_sea Vsp_sea];
 
 ng1_12k = [0 n1_12k(i112k:length(n1_12k)) vn_envelope(2,7:12) n6_12k(1:i612k) 0 0];
 vg1_12k = [Vsp_12k V_112k(i1sea:length(V_1sea)) vn_envelope(1,7:12) V_612k(1:i6sea) Vsn_12k Vsp_12k];
 
 plot(v_points_sea, n_points_sea,'b')
 hold on
 plot(vn_gust(1:6,1),vn_gust(1:6,2),'r--')
 hold on
 plot(vg1_sea,ng1_sea,'g')
 grid on
 xlabel('velocity (m/s)')
 ylabel('load factor (gs)')
 title(' vn diagram sea level')
 legend('normal','gust','envelope')
 
 figure (2)
 plot(v_points_12k, n_points_12k,'b')
 hold on
 plot(vn_gust(7:12,1),vn_gust(7:12,2),'r--')
 hold on
 plot(vg1_12k,ng1_12k,'g')
 grid on
 xlabel('velocity (m/s)')
 ylabel('load factor (gs)')
 title(' vn diagram 12k feet')
 legend('normal','gust','envelope')          