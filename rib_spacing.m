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
A_stringer = 9.8e-5; %[m^2]
% Ix = 3.1707e-3;     %[m^4]
% Iy = 1.381e-1;      %[m^4]
% Ixy = 1.9783e-3;     %[m^4]

% assuming stringer is a cirular rod
I = 7.6426e-10;
z = linspace(0,5.6,560);

finished = false;
k = 1;
i = 1;
L_t = 0;
while ~finished
    Pcr = A_stringer*abs(min_z(k))*1.5;
    L(i) = 2*sqrt((pi^2)*E*I/Pcr);
    L_t = sum(L);
    if L_t > 5.6
        finished = true;
    else
        k = find(abs(z - L_t) < 0.005);
        i = i + 1;
    
        if k > 560 || i > 20
            break
        end
    end
    
end

L_z = zeros(1,9);
for i =1:length(L)
    if i == 1
        L_z(i) = 0;
    else
        L_z(i) = sum(L(1:i - 1));
    end
end

Y = [1.3/2 1.3/2 -1.3/2 -1.3/2];
Z = [0 5.6 5.6 0];
plot(Z,Y,'b')
hold on
for i = 1:length(L_z)
    plot([L_z(i) L_z(i)],[-1.3/2 1.3/2],'r')
    hold on
end
legend('wing','ribs')
xlabel('spanwise direction (m)')
ylabel('chord direction (m)')

