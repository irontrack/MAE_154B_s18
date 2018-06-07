% shear buckling code 

% unfinished, when it is finished, it should load 
% all panels and compare all panels between each rib.
load('Sigma_ZZ.mat')
load('TopStringers.mat')
load('BotStringers.mat')
 X = TopStringers.posX;
 Y = TopStringers.posY;
 dX = diff(X);
 dY = diff(Y);
 B = sqrt(dX.^2 + dY.^2);
 
sigma_z = zeros(12,560);

 % cycle through each cell and find the maximum compressive stress
 % 
for i =1:12
    
    for j =1:560
        temp = SigmaZ{i,j};
        sigma_z(i,j) = min(min(temp));
        
    end
    
end
c = 1.346;
min_z = min(sigma_z);
% b = sqrt((0.9381 - 0.7478)^2 + (-0.0218 - -0.0305)^2) % length of longest panel

v = 0.33;
E = 73.1e9;
k = 11;
t = 0.001;
for i = 1:length(B)
    b = B(i);
    sigma_cr = (k*(pi^2)*E*t^2)/(12*(1 - v^2)*b^2)

    if abs(min_z) < sigma_cr
        fprintf('passed\n')
        fos = sigma_cr/abs(min(min_z));
        fprintf('fos: %4.2f\n',fos)
    else
        fprintf('not passed\n')
    end
end
