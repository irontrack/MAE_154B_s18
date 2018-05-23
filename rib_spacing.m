% BUCKLING LOAD AND RIB LENGTH PROGRAM

load('Sigma_ZZ.mat')

sigma_z = zeros(12,560);

 % cycle through each cell and find the maximum compressive stress
 % 
for i =1:12
    
    for j =1:560
        temp = SigmaZ{i,j};
        sigma_z(i,j) = min(min(temp));
        
    end
    
end

min_z = min(sigma_z);


% using ALUMINUM 2024 properties
E = 73.1e9;
A_stringer = 0.098; %[m^2]
% Ix = 3.1707e-3;     %[m^4]
% Iy = 1.381e-1;      %[m^4]
% Ixy = 1.9783e-3;     %[m^4]

% assuming stringer is a cirular rod
I = 7.6426e-04;
z = linspace(0,5.6,560);

finished = false;
k = 1;
i = 1;
while ~finished
    Pcr = A_stringer*abs(min_z(k))*1.5;
    Le(i) = sqrt((pi^2)*E*I/Pcr);
    L_t = sum(Le);
    k = find(abs(z - L_t) < 0.01);
    i = i + 1;
    
    if k > 560 || i > 20
        break
    end
    
    
end
