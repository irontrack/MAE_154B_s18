%The first part is the same as Airfoil_Structural_Setup.m. 
%Airfoil set-up (component dimensions assignment, airfoil contour graphing , 
%centroid, moments of inertia calcs). 
%This script utilizes the "struct" variable type (the "." notation) in
%MATLAB.
%Download the get_z function from CCLE homepage
close all;

clc;

chord=5; %chord length

Vx = 1; %test loads will be applied individually
%Vz = 1; My = 1;  

%define a few stringers (divide the airfoil into 4 segments: top, bottom,
%nose top, nose bottom)
numTopStringers = 5;
numBottomStringers = 5;
numNoseTopStringers = 5;
numNoseBottomStringers = 5;

%Web skin thickness
t_upper = 0.02/12;
t_lower = 0.02/12;
t_upper_front = 0.02/12;
t_lower_front = 0.02/12;
t_frontSpar = 0.04/12;
t_rearSpar = 0.04/12;

%Coordinate system definition:
%x-dir is from nose to tail, z-dir is normal to thw wing, y-dir is root to
%tip

%Positioning of structural elements
frontSpar = 0.2*chord;
backSpar = 0.7*chord;

% needs to be updated with our own geometry
% rename all sparcap struct usages
SparCaps.posX(1) = frontSpar;  %x location
sparCaps(2).posX = frontSpar;
sparCaps(3).posX = backSpar;
sparCaps(4).posX = backSpar;

sparCaps(1).posZ = chord*get_z(frontSpar/chord,1); %download the get_z function from CCLE
sparCaps(2).posZ = chord*get_z(frontSpar/chord,0); %z location
sparCaps(3).posZ = chord*get_z(backSpar/chord,1);
sparCaps(4).posZ = chord*get_z(backSpar/chord,0);

%Sizing of structural elements
sparCaps(1).areaO = .5; %area
sparCaps(2).areaO = .5;
sparCaps(3).areaO = .5;
sparCaps(4).areaO = .5;

%Gaps between adjacent stringers and spar caps
TopStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numTopStringers + 1);
BottomStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numBottomStringers + 1);
NoseTopStringerGap = (sparCaps(1).posX - 0)/(numNoseTopStringers + 1);
NoseBottomStringerGap = (sparCaps(1).posX - 0)/(numNoseBottomStringers + 1);


%Assume stringers spaced evenly along X axis betwen Spars in the following.
%Positioning and sizing of structural elements:

%top Stringers
for i=1:numTopStringers
    TopStringers(i).posX = sparCaps(1).posX + TopStringerGap*i;
    TopStringers(i).posZ = chord*get_z(TopStringers(i).posX/chord,1); %download the get_z function from CCLE
    TopStringers(i).areaO = 0.1; %area
end

%bottom Stringers
for i=1:numBottomStringers
    BottomStringers(i).posX = sparCaps(4).posX - BottomStringerGap*i;
    BottomStringers(i).posZ = chord*get_z(BottomStringers(i).posX/chord,0);
    BottomStringers(i).areaO = 0.1;

end

%nose bottom Stringers
for i=1:numNoseBottomStringers
    NoseBottomStringers(i).posX = sparCaps(2).posX - NoseBottomStringerGap*i;
    NoseBottomStringers(i).posZ = chord*get_z(NoseBottomStringers(i).posX/chord,0);
    NoseBottomStringers(i).areaO = 0.1;
end

%nose top Stringers
for i=1:numNoseTopStringers
    NoseTopStringers(i).posX = NoseTopStringerGap*i;
    NoseTopStringers(i).posZ = chord*get_z(NoseTopStringers(i).posX/chord,1);
    NoseTopStringers(i).areaO = 0.1;
end

%Compute the airfoil centroid.
%Calculate the contribution of stringers and sparcaps to airfoil centroid :
centroid.posX = sum([sparCaps.posX].*[sparCaps.areaO]) + ...
    sum([TopStringers.posX].*[TopStringers.areaO]) + ...
    sum([BottomStringers.posX].*[BottomStringers.areaO]) + ...
    sum([NoseTopStringers.posX].*[NoseTopStringers.areaO]) + ...
    sum([NoseBottomStringers.posX].*[NoseBottomStringers.areaO]);

centroid.posX = centroid.posX / ( sum([sparCaps.areaO]) + sum([TopStringers.areaO]) + ...
    sum([BottomStringers.areaO]) + sum([NoseTopStringers.areaO]) + sum([NoseBottomStringers.areaO]));

%%% ToDo Task: calculate the z-coordinate of the centroid.
% centroid.posZ = ...
%%%

%%% ToDo Task: 
% calculate the contribution of spars and skins to the centroid of the
% airfoil, and add it to the centroid calculated immediatley above. 
% Hint: do centroid.posX and centroid.posZ just like shown above.
%%%


%Plotting airfoil cross-section contour w/ centroid 
xChord = 0:.01:1*chord;
upperSurface = zeros(1,length(xChord));
lowerSurface = zeros(1,length(xChord));

for i=1:length(xChord)
    upperSurface(i) = chord*get_z(xChord(i)/chord,1);
    lowerSurface(i) = chord*get_z(xChord(i)/chord,0);
end

figure; hold on; axis equal; grid on;
plot(xChord,0*xChord,'-')
plot(xChord,upperSurface,'-')
plot(xChord,lowerSurface,'-')

plot([sparCaps(1).posX sparCaps(2).posX],[sparCaps(1).posZ sparCaps(2).posZ],'-')
plot([sparCaps(3).posX sparCaps(4).posX],[sparCaps(3).posZ sparCaps(4).posZ],'-')
plot([sparCaps.posX],[sparCaps.posZ],'o')
plot([TopStringers.posX],[TopStringers.posZ],'or')
plot([BottomStringers.posX],[BottomStringers.posZ],'or')
plot([NoseTopStringers.posX],[NoseTopStringers.posZ],'or')
plot([NoseBottomStringers.posX],[NoseBottomStringers.posZ],'or')
plot(centroid.posX,centroid.posZ,'rx') %centroid due to stringers, sparcaps, spars and skin
title('NACA 2412 cross-section w/ centroid marked as "x"')


% Summing contributions for inertia terms
Ix = 0; Iz = 0; Ixz = 0; %moments of inertia

for i=1:4 %spar caps
    Ix = Ix + sparCaps(i).areaO*(sparCaps(i).posZ-centroid.posZ)^2;
    Iz = Iz + sparCaps(i).areaO*(sparCaps(i).posX-centroid.posX)^2;
    Ixz = Ixz + sparCaps(i).areaO*(sparCaps(i).posX-centroid.posX)*(sparCaps(i).posZ-centroid.posZ);
end

for i=1:numTopStringers %top stringers
    Ix = Ix + TopStringers(i).areaO*(TopStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + TopStringers(i).areaO*(TopStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + TopStringers(i).areaO*(TopStringers(i).posX-centroid.posX)*(TopStringers(i).posZ-centroid.posZ);
end
for i=1:numBottomStringers %bottom stringers
    Ix = Ix + BottomStringers(i).areaO*(BottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + BottomStringers(i).areaO*(BottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + BottomStringers(i).areaO*(BottomStringers(i).posX-centroid.posX)*(BottomStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseTopStringers %nose top stringers
    Ix = Ix + NoseTopStringers(i).areaO*(NoseTopStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + NoseTopStringers(i).areaO*(NoseTopStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + NoseTopStringers(i).areaO*(NoseTopStringers(i).posX-centroid.posX)*(NoseTopStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseBottomStringers %nose bottom stringers
    Ix = Ix + NoseBottomStringers(i).areaO*(NoseBottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + NoseBottomStringers(i).areaO*(NoseBottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + NoseBottomStringers(i).areaO*(NoseBottomStringers(i).posX-centroid.posX)*(NoseBottomStringers(i).posZ-centroid.posZ);
end

%%% ToDo Task: 
% calculate the contribution of spars and skins to the moments of inertia of the airfoil,
% and add it to the inertia terms calculated immediatley above. 
%%%


%stacking all the stringers and spar caps
stringersTotal=[NoseTopStringers,sparCaps(1),TopStringers,sparCaps(3),NoseBottomStringers,BottomStringers,sparCaps(2),sparCaps(4)];



%% Skin panels creation and boom area adjustment (structural idealization)
%%Download the get_ds.m and get_Bs.m functions from CCLE homepage

%Airfoil cell 1 (right cell)

%upper webs
numStringers = numTopStringers;
stringerGap = TopStringerGap;
webThickness = t_upper;
tempStringers = TopStringers;

for i=1:(numStringers+1)
    %web = the skin panel between 2 adjacent stringers or between 1 stringer and its adjacent sparcap. 
    web(i).xStart = sparCaps(1).posX + stringerGap*(i-1); %assigning x abd z coordinates to the ends of webs
    web(i).xEnd = sparCaps(1).posX + stringerGap*(i); 
    web(i).thickness = webThickness;
    web(i).zStart = chord*get_z(web(i).xStart/chord,1);
    web(i).zEnd = chord*get_z(web(i).xEnd/chord,1);
    web(i).ds = chord*get_ds(web(i).xStart/chord,web(i).xEnd/chord,1); %web length
                                    %download the get_ds function from CCLE
    
   %structural idealization
   if i==1 %first web
        B=get_Bs(sparCaps(1).posZ,TopStringers(i).posZ,t_upper,centroid.posZ,web(i).ds); %download the get_Bs function from CCLE
        sparCaps(1).area=sparCaps(1).areaO+B(1); %transferring the skin thickness to areas of adjacent stringers or sparcaps
        TopStringers(i).area=TopStringers(i).areaO+B(2);
    elseif i==numTopStringers+1 %last web
        B=get_Bs(TopStringers(i-1).posZ,sparCaps(3).posZ,t_upper,centroid.posZ,web(i).ds); 
        sparCaps(3).area=sparCaps(3).areaO+B(2);
        TopStringers(i-1).area=TopStringers(i-1).areaO+B(1);
   else %anything in between
        B=get_Bs(TopStringers(i-1).posZ,TopStringers(i).posZ,t_upper,centroid.posZ,web(i).ds);
        TopStringers(i-1).area=TopStringers(i-1).areaO+B(1);
        TopStringers(i).area=TopStringers(i).areaO+B(2);
    end
end
webTop = web; 
web = [];

%rear spar
i=1;
web(i).xStart = sparCaps(3).posX; %sparcap3 is the top right spar cap
web(i).xEnd = sparCaps(4).posX; %sparcap 4 is the bottom right spar cap
web(i).thickness = t_rearSpar;
web(i).zStart = sparCaps(3).posZ;
web(i).zEnd = sparCaps(4).posZ;
web(i).ds = abs(sparCaps(3).posZ - sparCaps(4).posZ); 

B=get_Bs(sparCaps(3).posZ,sparCaps(4).posZ,t_rearSpar,centroid.posZ, sparCaps(3).posZ-sparCaps(4).posZ);   
sparCaps(3).area=sparCaps(3).areaO+B(1);
sparCaps(4).area=sparCaps(4).areaO+B(2);

webRearSpar = web;
web = [];


%lower webs
numStringers = numBottomStringers;
stringerGap = BottomStringerGap;
webThickness = t_lower;
tempStringers = BottomStringers;

%%% ToDo Task:
for i=1:(numStringers+1)
    %web(i).xStart = ...
    %web(i).xEnd = ...
    %web(i).thickness = ...
    %web(i).zStart = ...
    %web(i).zEnd = ...
    %web(i).ds = ...

    if i==1  %first web
        %B= ...
        %sparCaps(2).area=...
        %BottomStringers(i).area=...
    elseif i==numBottomStringers+1  %last web
        %B=...
        %sparCaps(4).area=...
        %BottomStringers(i-1).area=...
    else  %anything in between
        %B=...
        %BottomStringers(i-1).area=...
        %BottomStringers(i).area=...
    end
end
%%%
webBottom = web;
web = [];


%front Spar
i=1;
web(i).xStart = sparCaps(2).posX; %sparcap2 is the bottom left spar cap
web(i).xEnd = sparCaps(1).posX; %sparcap1 is the top left spar cap
web(i).thickness = t_frontSpar;
web(i).zStart = sparCaps(2).posZ;
web(i).zEnd = sparCaps(1).posZ;
web(i).ds = abs(sparCaps(2).posZ - sparCaps(1).posZ);
%%% ToDo Task:
%B=...
%sparCaps(1).area=...
%sparCaps(2).area=...
%%%
webFrontSpar = web;
web = [];


%% Airfoil cell 2 (left cell)

%lower nose webs
numStringers = numNoseBottomStringers;
stringerGap = NoseBottomStringerGap;
webThickness = t_lower_front;
tempStringers = NoseBottomStringers;

%%% ToDo Task:
for i=1:(numStringers+1)
    %web(i).xStart = ...
    %web(i).xEnd = ...
    %web(i).thickness = ...
    %web(i).zStart = ...
    %web(i).zEnd = ...
    %web(i).ds = ...
  
    if i==numNoseBottomStringers+1 %last web
        %B=...
        %sparCaps(2).area=...
        %NoseBottomStringers(i-1).area=...
    elseif i==1 %first web
        NoseBottomStringers(i).area=NoseBottomStringers(i).areaO+t_lower_front*web(i).ds/3; 
    else %anything in between
        %B=...
        %NoseBottomStringers(i-1).area=...
        %NoseBottomStringers(i).area=...
    end
end
%%%
webLowerNose = web;
web = [];

%upper nose webs
numStringers = numNoseTopStringers;
stringerGap = NoseTopStringerGap;
webThickness = t_upper_front;
tempStringers = NoseTopStringers;

%%% ToDo Task:
for i=1:(numStringers+1)
    %web(i).xStart = ...
    %web(i).xEnd = ...
    %web(i).thickness = ...
    %web(i).zStart = ...
    %web(i).zEnd = ...
    %web(i).ds = ...
  
   if i==numNoseTopStringers+1
        %B=...
        %sparCaps(1).area=...
        %NoseTopStringers(i-1).area=...
    elseif i==1
        NoseTopStringers(i).area=NoseTopStringers(i).areaO+t_upper_front*web(i).ds/3; 
    else
        %B=...
        %NoseTopStringers(i-1).area=...
        %NoseTopStringers(i).area=...
   end
end
%%%
webUpperNose = web;
web = [];


%front Spar
i=1;
web(i).xStart = sparCaps(1).posX;
web(i).xEnd = sparCaps(2).posX;
web(i).thickness = t_frontSpar;
web(i).zStart = sparCaps(1).posZ;
web(i).zEnd = sparCaps(2).posZ;
web(i).ds = abs(sparCaps(1).posZ - sparCaps(2).posZ);

webFrontSparCell2 = web;
web = [];

%% recalculate centroid
centroid.posX = sum([sparCaps.posX].*[sparCaps.area]) + ...
    sum([TopStringers.posX].*[TopStringers.area]) + ...
    sum([BottomStringers.posX].*[BottomStringers.area]) + ...
    sum([NoseTopStringers.posX].*[NoseTopStringers.area]) + ...
    sum([NoseBottomStringers.posX].*[NoseBottomStringers.area]);

centroid.posX = centroid.posX / ( sum([sparCaps.area]) + sum([TopStringers.area]) + ...
    sum([BottomStringers.area]) + sum([NoseTopStringers.area]) + sum([NoseBottomStringers.area]));

%%% ToDo Task:
% centroid.posZ = ...
%%%   


%% resumming contributions for inertia terms
Ix = 0; Iz = 0; Ixz = 0;

for i=1:4 %spar caps
    Ix = Ix + sparCaps(i).area*(sparCaps(i).posZ-centroid.posZ)^2;
    Iz = Iz + sparCaps(i).area*(sparCaps(i).posX-centroid.posX)^2;
    Ixz = Ixz + sparCaps(i).area*(sparCaps(i).posX-centroid.posX)*(sparCaps(i).posZ-centroid.posZ);
end


for i=1:numTopStringers %top stringers
    Ix = Ix + TopStringers(i).area*(TopStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + TopStringers(i).area*(TopStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + TopStringers(i).area*(TopStringers(i).posX-centroid.posX)*(TopStringers(i).posZ-centroid.posZ);
end
for i=1:numBottomStringers %bottom stringers
    Ix = Ix + BottomStringers(i).area*(BottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + BottomStringers(i).area*(BottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + BottomStringers(i).area*(BottomStringers(i).posX-centroid.posX)*(BottomStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseTopStringers %nose top stringers
    Ix = Ix + NoseTopStringers(i).area*(NoseTopStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + NoseTopStringers(i).area*(NoseTopStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + NoseTopStringers(i).area*(NoseTopStringers(i).posX-centroid.posX)*(NoseTopStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseBottomStringers %nose bottom stringers
    Ix = Ix + NoseBottomStringers(i).area*(NoseBottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + NoseBottomStringers(i).area*(NoseBottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + NoseBottomStringers(i).area*(NoseBottomStringers(i).posX-centroid.posX)*(NoseBottomStringers(i).posZ-centroid.posZ);
end


%stacking all the stringers and spar caps
stringersTotal=[NoseTopStringers,sparCaps(1),TopStringers,sparCaps(3),NoseBottomStringers,BottomStringers,sparCaps(2),sparCaps(4)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The code below this point are new to you guys (for week 5 onward). You
%should have finished everything above this point. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shear flow calculation ---- please read Shear_Flow_MATLAB_Guide.pdf posted on CCLE
%Download all the functions from the homepage on CCLE (get_int.m, get_dp.m,
%etc.)

%Cell 1 (right cell)

%upper webs
numStringers = numTopStringers;
stringerGap = TopStringerGap;
webThickness = t_upper;
tempStringers = TopStringers;

for i=1:(numStringers+1)
    web(i).xStart = sparCaps(1).posX + stringerGap*(i-1);
    web(i).xEnd = sparCaps(1).posX + stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = chord*get_z(web(i).xStart/chord,1);
    web(i).zEnd = chord*get_z(web(i).xEnd/chord,1);
    if i==1 %the left-most spar cap
        web(i).dp_area = sparCaps(1).area;
        web(i).dP_X = 0; %set delta_p to be zero in the first spar cap
        web(i).dP_Z = 0;
        web(i).qPrime_X = 0;
        web(i).qPrime_Z = 0;
    else %all stringers
        web(i).dp_area = tempStringers(i-1).area; %stringer area
        dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ; %web location
        web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);  %delta_P due to just Vx shear
        web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);  %....just Vz shear
        web(i).qPrime_X = web(i-1).qPrime_X - web(i).dP_X; %accumulating delta_p terms across all webs 
        web(i).qPrime_Z = web(i-1).qPrime_Z - web(i).dP_Z; %accumulating.....
    end
    tempInt = chord*get_int(web(i).xStart/chord,web(i).xEnd/chord,1);  %integral of airfoil profile function
    triangle1 = abs( (web(i).xStart - sparCaps(1).posX)*web(i).zStart/2); %geometric trick (pls try to understand this)
    triangle2 = abs((web(i).xEnd - sparCaps(1).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2; %finding area swept out by webs
    web(i).ds = chord*get_ds(web(i).xStart/chord,web(i).xEnd/chord,1); %web length
    web(i).dS_over_t = web(i).ds / web(i).thickness; %web length divided by web thickness
    
    %Preparing terms for the twist equations and the torque equations
    web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t; %twist terms
    web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
    web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X; %torque terms
    web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
    web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart); %dimensionalize the terms by length (q'*dx)
    web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
    web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
    web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);
    web(i).radCurv = get_curve(web(i)); %radius of curvature of web
end
webTop = web;
web = [];

%rear spar
i=1;
web(i).xStart = sparCaps(3).posX;
web(i).xEnd = sparCaps(4).posX;
web(i).thickness = t_rearSpar;
web(i).zStart = sparCaps(3).posZ;
web(i).zEnd = sparCaps(4).posZ;
web(i).dp_area = sparCaps(3).area;
dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime_X = webTop(numTopStringers+1).qPrime_X - web(i).dP_X;
web(i).qPrime_Z = webTop(numTopStringers+1).qPrime_Z - web(i).dP_Z;

web(i).Area = (sparCaps(3).posX-sparCaps(1).posX)*sparCaps(3).posZ/2 + ...
    abs((sparCaps(3).posX-sparCaps(1).posX)*sparCaps(4).posZ/2);
web(i).ds = abs(sparCaps(3).posZ - sparCaps(4).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;

web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);
web(i).radCurv=1;
webRearSpar = web;
web = [];


%lower webs
numStringers = numBottomStringers;
stringerGap = BottomStringerGap;
webThickness = t_lower;
tempStringers = BottomStringers;

for i=1:(numStringers+1)
   %%% 
   %ToDo Task
   %%% 
end
webBottom = web;
web = [];


%front Spar
i=1;
   %%% 
   %ToDo Task
   %%% 
web(i).radCurv=1;
webFrontSpar = web;
web = [];




%% Cell 2 (left cell)

%lower nose webs
numStringers = numNoseBottomStringers;
stringerGap = NoseBottomStringerGap;
webThickness = t_lower_front;
tempStringers = NoseBottomStringers;

for i=1:(numStringers+1)
    
   %%% 
   %ToDo Task
   %%% 

end
webLowerNose = web;
web = [];


%upper nose webs
numStringers = numNoseTopStringers;
stringerGap = NoseTopStringerGap;
webThickness = t_upper_front;
tempStringers = NoseTopStringers;

for i=1:(numStringers+1)
    
   %%% 
   %ToDo Task
   %%% 
end
webUpperNose = web;
web = [];


%front Spar
i=1;
   %%% 
   %ToDo Task
   %%% 
web(i).radCurv=1;
webFrontSparCell2 = web;
web = [];


%Code validaiton:
%check that q'*dx sums up to Vx

Fx = sum([webTop.qp_dx_X])+webRearSpar.qp_dx_X+ sum([webBottom.qp_dx_X])+webFrontSpar.qp_dx_X;  %cell 1
Fx = Fx + sum([webLowerNose.qp_dx_X])+ sum([webUpperNose.qp_dx_X]);  %cell 2
Fx

% Fz = sum([webTop.qp_dz_X])+webRearSpar.qp_dz_X+ sum([webBottom.qp_dz_X])+webFrontSpar.qp_dz_X;  %cell 1
% Fz = Fz + sum([webLowerNose.qp_dz_X])+ sum([webUpperNose.qp_dz_X]);  %cell 2
% Fz

%% Shear flow solution procedure----pure bending + pure torsion
%With all parameters set up, we can now start solving...

% Sum up the ds/t and q*ds/t to solve 2 twist equations, 2 unknowns (pure
% bending) 

% Method: convert equations into matrix form: [A]*[q1s q2s] = B

A11 = sum([webTop.dS_over_t])+webRearSpar.dS_over_t+ sum([webBottom.dS_over_t])+webFrontSpar.dS_over_t;
A22 = sum([webLowerNose.dS_over_t])+ sum([webUpperNose.dS_over_t])+webFrontSparCell2.dS_over_t;
A12 = -webFrontSpar.dS_over_t;
A21 = -webFrontSparCell2.dS_over_t;

B1_X = sum([webTop.q_dS_over_t_X])+webRearSpar.q_dS_over_t_X+ sum([webBottom.q_dS_over_t_X])+webFrontSpar.q_dS_over_t_X;
B2_X = sum([webLowerNose.q_dS_over_t_X])+ sum([webUpperNose.q_dS_over_t_X])+webFrontSparCell2.q_dS_over_t_X;
B1_Z = sum([webTop.q_dS_over_t_Z])+webRearSpar.q_dS_over_t_Z+ sum([webBottom.q_dS_over_t_Z])+webFrontSpar.q_dS_over_t_Z;
B2_Z = sum([webLowerNose.q_dS_over_t_Z])+ sum([webUpperNose.q_dS_over_t_Z])+webFrontSparCell2.q_dS_over_t_Z;

Amat = [A11 A12; A21 A22];
Bmat_X = -[B1_X;B2_X];
Bmat_Z = -[B1_Z;B2_Z];

qs_X = inv(Amat)*Bmat_X; %pure bending solution due to Vx
qs_Z = inv(Amat)*Bmat_Z; %pure bending solution due to Vz


%Finding the shear center
sum_2_a_q_X = sum([webTop.two_A_qprime_X])+webRearSpar.two_A_qprime_X+ sum([webBottom.two_A_qprime_X]);  %cell 1 qprimes
sum_2_a_q_X = sum_2_a_q_X + sum([webLowerNose.two_A_qprime_X])+ sum([webUpperNose.two_A_qprime_X]);   %cell 2 qprimes
sum_2_a_q_X = sum_2_a_q_X +  qs_X(1)*(sum([webTop.Area])+webRearSpar.Area+ sum([webBottom.Area]));
sum_2_a_q_X = sum_2_a_q_X +  qs_X(2)*(sum([webLowerNose.Area])+ sum([webUpperNose.Area]));
%%% ToDo Task:
%sum_2_a_q_Z = ...             %cell 1 qprimes
%sum_2_a_q_Z = ...             %cell 2 qprimes
%sum_2_a_q_Z = ...
%sum_2_a_q_Z = ...
%%%

%shear center
sc.posX =  sum_2_a_q_Z / Vz + frontSpar;
sc.posZ =  sum_2_a_q_X / Vx;


% Now consider the torque representing shifting the load from the quarter
% chord to the SC  (need to check signs on these moments)
 
torque_Z = Vz*(sc.posX - 0.25);
torque_X = Vx*sc.posZ;

Area1 = sum([webTop.Area]) + webRearSpar.Area + sum([webBottom.Area]);
%check area
Area1_check = chord*get_int(frontSpar/chord,backSpar/chord,1) + chord*get_int(frontSpar/chord,backSpar/chord,0);

Area2 = sum([webLowerNose.Area]) + sum([webUpperNose.Area]);
Area2_check = chord*get_int(0,frontSpar/chord,1) + chord*get_int(0,frontSpar/chord,0);


%For twist equation (pure torsion):

q1t_over_q2t = (A22/Area2 + webFrontSpar.dS_over_t/Area1)/(A11/Area1 + webFrontSpar.dS_over_t/Area2);

q2t = torque_X/(2*Area1*q1t_over_q2t + 2*Area2);
q1t = q2t*q1t_over_q2t;
qt_X = [q1t;q2t]; %pure torsion solution due to Vx

%%% ToDo Task:
%q2t = 
%q1t = 
%qt_Z =           %pure torsion solution due to Vz
%%%

q2t = My/(2*Area1*q1t_over_q2t + 2*Area2);
q1t = q2t*q1t_over_q2t;
qt_My = [q1t;q2t]; %pure torsion solution due to My


%Adding up qs, qt, and qprime terms to get total shear flow along webs
for i=1:numBottomStringers+1
webBottom(i).Vx_tau=(webBottom(i).qPrime_X+qt_X(1)+qs_X(1))/t_lower;
webBottom(i).Vz_tau=(webBottom(i).qPrime_Z+qt_Z(1)+qs_Z(1))/t_lower;
webBottom(i).My_tau = qt_My(1)/t_lower;
end

for i=1:numTopStringers+1
%%% ToDo Task %%%
end

for i=1:numNoseBottomStringers+1
%%% ToDo Task %%%
end

for i=1:numNoseTopStringers+1
webUpperNose(i).Vx_tau=(webUpperNose(i).qPrime_X+qt_X(2)+qs_X(2))/t_upper_front;
webUpperNose(i).Vz_tau=(webUpperNose(i).qPrime_Z+qt_Z(2)+qs_Z(2))/t_upper_front;
webUpperNose(i).My_tau = qt_My(2)/t_upper_front;
end

webFrontSpar.Vx_tau=(webFrontSpar.qPrime_X+qt_X(1)+qs_X(1)-qt_X(2)-qs_X(2))/t_frontSpar;
webFrontSpar.Vz_tau=(webFrontSpar.qPrime_Z+qt_Z(1)+qs_Z(1)-qt_Z(2)-qs_Z(2))/t_frontSpar;
webFrontSpar.My_tau=qt_My(1)/t_frontSpar;

webFrontSparCell2.Vx_tau=(webFrontSparCell2.qPrime_X+qt_X(2)+qs_X(2)-qt_X(1)-qs_X(1))/t_frontSpar;
webFrontSparCell2.Vz_tau=(webFrontSparCell2.qPrime_Z+qt_Z(2)+qs_Z(2)-qt_Z(1)-qs_Z(1))/t_frontSpar;
webFrontSparCell2.My_tau=qt_My(1)/t_frontSpar;

%%% ToDo Task
%webRearSpar.Vx_tau=...
%webRearSpar.Vz_tau=...
%webRearSpar.My_tau=...
%%%

% stacking all the webs
webTotal=[webUpperNose,webTop,webLowerNose,webBottom,webFrontSpar,webRearSpar];

%stacking all the stringers
stringersTotal=[NoseTopStringers,sparCaps(1),TopStringers,sparCaps(3),NoseBottomStringers,BottomStringers,sparCaps(2),sparCaps(4)];


%plotting airfoil cross-section
xChord = 0:.01:1*chord;
upperSurface = zeros(1,length(xChord));
lowerSurface = zeros(1,length(xChord));

for i=1:length(xChord)
    upperSurface(i) = chord*get_z(xChord(i)/chord,1);
    lowerSurface(i) = chord*get_z(xChord(i)/chord,0);
end

figure; hold on; axis equal; grid on;
plot(xChord,0*xChord,'-')
plot(xChord,upperSurface,'-')
plot(xChord,lowerSurface,'-')

plot([sparCaps(1).posX sparCaps(2).posX],[sparCaps(1).posZ sparCaps(2).posZ],'-')
plot([sparCaps(3).posX sparCaps(4).posX],[sparCaps(3).posZ sparCaps(4).posZ],'-')
plot([sparCaps.posX],[sparCaps.posZ],'o')
plot([TopStringers.posX],[TopStringers.posZ],'or')
plot([BottomStringers.posX],[BottomStringers.posZ],'or')
plot([NoseTopStringers.posX],[NoseTopStringers.posZ],'or')
plot([NoseBottomStringers.posX],[NoseBottomStringers.posZ],'or')
plot(centroid.posX,centroid.posZ,'rx')
plot(sc.posX,sc.posZ,'gx')