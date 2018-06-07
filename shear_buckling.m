% shear buckling code 

% unfinished, when it is finished, it should load 
% all panels and compare all panels between each rib.

% b = sqrt((0.9381 - 0.7478)^2 + (-0.0218 - -0.0305)^2) % length of longest panel
b = 0.1;
v = 0.33;
E = 73.1e9;
k = 11;
t = 0.001;

sigma_cr = (k*(pi^2)*E*t^2)/(12*(1 - v^2)*b^2);