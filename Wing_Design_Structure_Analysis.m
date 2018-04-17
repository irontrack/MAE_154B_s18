%% MAE 154B Wing Design Structure Analysis
close all;
clear all;
clc;

%% Define Nodes
% Import Shape of NACA2412
TopNodes = importdata('NACA2412_Top.txt');
BottomNodes = importdata('NACA2412_Bottom.txt');
% Extract Node Points
Nx=TopNodes(:,1);
Nytop=TopNodes(:,2);
Nybot=BottomNodes(:,2);
% Scale points to desire size
scalefactor=1.346;
nx=scalefactor*Nx;
nytop=scalefactor*Nytop;
nybot=scalefactor*Nybot;

%% Create Elements
% All Nodes for airfoil are in the inner side of the elements
skint=0.003;  % Skin thickness

% Create Airfoil Elements
Topel=struct('posX',nx,'posY',nytop,'Area',zeros(length(nx),1),...
             'cx',zeros(length(nx),1),'cy',zeros(length(nx),1),...
             'Ixx',zeros(length(nx),1),'Iyy',zeros(length(nx),1),...
             'Ixy',zeros(length(nx),1),...
             'xcoord',zeros(length(nx),5),'ycoord',zeros(length(nx),5));
Botel=struct('posX',nx,'posY',nybot,'Area',zeros(length(nx),1),...
             'cx',zeros(length(nx),1),'cy',zeros(length(nx),1),...
             'Ixx',zeros(length(nx),1),'Iyy',zeros(length(nx),1),...
             'Ixy',zeros(length(nx),1),...
             'xcoord',zeros(length(nx),5),'ycoord',zeros(length(nx),5));

% Calculate Aera and Centroid of Each Element, Plot Results
figure()
hold on
XiAiairf=0;
YiAiairf=0;
for i=1:length(nx)-1
   % Top elements
   Topel.xcoord(i,:)=[Topel.posX(i) Topel.posX(i)...
           Topel.posX(i+1) Topel.posX(i+1) Topel.posX(i)];
   Topel.ycoord(i,:)=[Topel.posY(i) Topel.posY(i)+skint ...
           Topel.posY(i+1)+skint Topel.posY(i+1) Topel.posY(i)];
   % Calculate Area
   Topel.Area(i)=skint*sqrt((Topel.posX(i+1)-Topel.posX(i))^2+...
                 (Topel.posY(i+1)-Topel.posY(i))^2);
   % Calculate Element Centroid
   Topel.cx(i)=0.5*((Topel.xcoord(i,1)+Topel.xcoord(i,2))/2+...
                    (Topel.xcoord(i,3)+Topel.xcoord(i,4))/2);
   Topel.cy(i)=0.5*((Topel.ycoord(i,1)+Topel.ycoord(i,2))/2+...
                    (Topel.ycoord(i,3)+Topel.ycoord(i,4))/2);
   % Plot    
   plot(Topel.xcoord(i,:),Topel.ycoord(i,:))
   plot(Topel.cx(i), Topel.cy(i),'r*')
   
   % Bottom elements
   Botel.xcoord(i,:)=[Botel.posX(i) Botel.posX(i)...
           Botel.posX(i+1) Botel.posX(i+1) Botel.posX(i)];
   Botel.ycoord(i,:)=[Botel.posY(i)-skint Botel.posY(i) ...
           Botel.posY(i+1) Botel.posY(i+1)-skint Botel.posY(i)-skint];
   % Calculate Area
   Botel.Area(i)=skint*sqrt((Botel.posX(i+1)-Botel.posX(i))^2+...
                 (Botel.posY(i+1)-Botel.posY(i))^2);
   % Calculate Element Centroid
   Botel.cx(i)=0.5*((Botel.xcoord(i,1)+Botel.xcoord(i,2))/2+...
                    (Botel.xcoord(i,3)+Botel.xcoord(i,4))/2);
   Botel.cy(i)=0.5*((Botel.ycoord(i,1)+Botel.ycoord(i,2))/2+...
                    (Botel.ycoord(i,3)+Botel.ycoord(i,4))/2);
   % Plot 
   plot(Botel.xcoord(i,:),Botel.ycoord(i,:))
   plot(Botel.cx(i), Botel.cy(i),'r*')
   
   XiAiairf=XiAiairf+Topel.cx(i)*Topel.Area(i)+Botel.cx(i)*Botel.Area(i);
   YiAiairf=YiAiairf+Topel.cy(i)*Topel.Area(i)+Botel.cy(i)*Botel.Area(i);
end
% Centroid of Airfoil
AirfoilArea=sum(Topel.Area(:,1))+sum(Botel.Area(:,1));
Cxaf=XiAiairf/AirfoilArea;
Cyaf=YiAiairf/AirfoilArea;

% Create Spars Elements
spart=0.003; % Spar thickness
sparcapt=0.002; % Spar cap size
SparIndex=[8 13];
Spars=struct('posX',zeros(length(SparIndex),1),...
             'posY',zeros(length(SparIndex),1),...
             'Area',zeros(length(SparIndex),1),...
             'cx',zeros(length(SparIndex),1), ...
             'cy',zeros(length(SparIndex),1),...
             'Ixx',zeros(length(SparIndex),1),...
             'Iyy',zeros(length(SparIndex),1),...
             'Ixy',zeros(length(SparIndex),1),...
             'xcoord',zeros(length(SparIndex),5),...
             'ycoord',zeros(length(SparIndex),5));

XiAiSpar=0;
YiAiSpar=0;
for i=1:length(SparIndex)
   ind=SparIndex(i);
   Spars.xcoord(i,:)=[nx(ind)-spart/2 nx(ind)-spart/2 ...
                      nx(ind)+spart/2 nx(ind)+spart/2 ...
                      nx(ind)-spart/2];
   Spars.ycoord(i,:)=[nybot(ind) nytop(ind) ...
                      nytop(ind) nybot(ind) nybot(ind)];
   % Calculate Area
   Spars.Area(i)=spart*(nytop(ind)-nybot(ind));
   % Calculate Element Centroid
   Spars.cx(i)=0.5*((Spars.xcoord(i,1)+Spars.xcoord(i,2))/2+...
                    (Spars.xcoord(i,3)+Spars.xcoord(i,4))/2);
   Spars.cy(i)=0.5*((Spars.ycoord(i,1)+Spars.ycoord(i,2))/2+...
                    (Spars.ycoord(i,3)+Spars.ycoord(i,4))/2);
   % Plot    
   plot(Spars.xcoord(i,:),Spars.ycoord(i,:))
   plot(Spars.cx(i), Spars.cy(i),'r*')
   
   XiAiSpar=XiAiSpar+Spars.cx(i)*Spars.Area(i);
   YiAiSpar=YiAiSpar+Spars.cy(i)*Spars.Area(i);
end
% Centroid of Spars
SparArea=sum(Spars.Area(:,1));
Cxspar=XiAiSpar/SparArea;
Cyspar=YiAiSpar/SparArea;

% Calculate Total Centroid
Cx=(XiAiairf+XiAiSpar)/(AirfoilArea+SparArea);
Cy=(YiAiairf+YiAiSpar)/(AirfoilArea+SparArea);

plot(Cx,Cy,'*')
xlim([-0.05,1.4])
ylim([-0.15,0.15])
xlabel('Length (m)')
ylabel('Length (m)')
grid on
hold off

%% Calculate Inertia



