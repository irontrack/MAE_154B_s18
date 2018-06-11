% BUCKLING LOAD AND RIB LENGTH PROGRAM
min_sigma_z = zeros(12,560);
% Cycle through each cell and find the maximum compressive stress
for LC =1:12
    for zi =1:560
        min_sigma_z(LC,zi) = min(min(SigmaZ{LC,zi}));     
    end
end
min_z = min(min_sigma_z);

% Assuming stringer is a z-shape
c = L_stringer - t_stringer;
I_Stringer=(L_stringer*(H_stringer^3)-c*(H_stringer-2*t_stringer)^3)/12;

finished = false;
k = 1;
i = 1;
L_t = 0;
while ~finished
    Pcr = StringerArea1*abs(min_z(k))*1.5;
    L(i) = 2*sqrt((pi^2)*E*I_Stringer/Pcr);
    L_t = sum(L);
    if L_t > 5.6
        finished = true;
    else
        k = min(find(abs(z - L_t) < 0.006));
        i = i + 1;
    
        if k > 560 || i > 30
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
figure()
plot(Z,Y,'b')
hold on
for i = 1:length(L_z)
    plot([L_z(i) L_z(i)],[-1.3/2 1.3/2],'r')
    hold on
end
legend('wing','ribs')
xlabel('spanwise direction (m)')
ylabel('chord direction (m)')