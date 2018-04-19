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
xscalefactor=1.346;
topscalefactor=1.325;
botscalefactor=1.38;
nx=xscalefactor*Nx;
nytop=topscalefactor*Nytop;
nybot=botscalefactor*Nybot;

%% Create Elements
% All Nodes for airfoil are in the inner side of the elements
skint=0.003;  % Skin thickness

% Create Airfoil Elements
Topel=struct('posX',nx,'posY',nytop,...
             'Length',zeros(length(nx),1),'Area',zeros(length(nx),1),...
             'cx',zeros(length(nx),1),'cy',zeros(length(nx),1),...
             'Ixx',zeros(length(nx),1),'Iyy',zeros(length(nx),1),...
             'Ixy',zeros(length(nx),1),...
             'xcoord',zeros(length(nx),5),'ycoord',zeros(length(nx),5));
Botel=struct('posX',nx,'posY',nybot,...
             'Length',zeros(length(nx),1),'Area',zeros(length(nx),1),...
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
   Topel.Length(i)=sqrt((Topel.posX(i+1)-Topel.posX(i))^2+...
                 (Topel.posY(i+1)-Topel.posY(i))^2);
   Topel.Area(i)=skint*Topel.Length(i);
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
   Botel.Length(i)=sqrt((Botel.posX(i+1)-Botel.posX(i))^2+...
                 (Botel.posY(i+1)-Botel.posY(i))^2);
   Botel.Area(i)=skint*Botel.Length(i);
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
             'Length',zeros(length(SparIndex),1),...
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
   Spars.Length(i)=(nytop(ind)-nybot(ind));
   Spars.Area(i)=spart*Spars.Length(i);
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
ylim([-0.5,0.5])
xlabel('Length (m)')
ylabel('Length (m)')
grid on
hold off

disp('The centroid is at')
disp([Cx,Cy])
%% Calculate Inertia

% For Airfoil Elements
for i=1:length(nx)-1
    %Top Element
    betatop=atan2(Topel.posY(i+1)-Topel.posY(i),...
                  Topel.posX(i+1)-Topel.posX(i));
    Topel.Ixx(i)=Topel.Length(i)^3*skint*(sin(betatop))^2/12+...
                 Topel.Area(i)*(Cy-Topel.cy(i))^2;
    Topel.Iyy(i)=Topel.Length(i)^3*skint*(cos(betatop))^2/12+...
                 Topel.Area(i)*(Cx-Topel.cx(i))^2;
    Topel.Ixy(i)=Topel.Length(i)^3*skint*sin(2*betatop)/24+...
                 Topel.Area(i)*(Cx-Topel.cx(i))*(Cy-Topel.cy(i));
    % Bottom Element
    betabot=atan2(Botel.posY(i+1)-Botel.posY(i),...
                  Botel.posX(i+1)-Botel.posX(i));
    Botel.Ixx(i)=Botel.Length(i)^3*skint*(sin(betabot))^2/12+...
                 Botel.Area(i)*(Cy-Botel.cy(i))^2;
    Botel.Iyy(i)=Botel.Length(i)^3*skint*(cos(betabot))^2/12+...
                 Botel.Area(i)*(Cx-Botel.cx(i))^2;
    Botel.Ixy(i)=Botel.Length(i)^3*skint*sin(2*betabot)/24+...
                 Botel.Area(i)*(Cx-Botel.cx(i))*(Cy-Botel.cy(i));
end

% Spars Elements
for i=1:length(SparIndex)
    Spars.Ixx(i)=spart*Spars.Length(i)^3/12+...
                 Spars.Area(i)*(Cy-Spars.cy(i))^2;
    Spars.Iyy(i)=spart^3*Spars.Length(i)/12+...
                 Spars.Area(i)*(Cx-Spars.cx(i))^2;
    Spars.Ixy(i)=Spars.Area(i)*(Cx-Spars.cx(i))*(Cy-Spars.cy(i));
end

Ixx=sum(Topel.Ixx(:))+sum(Botel.Ixx(:))+sum(Spars.Ixx(:));
Iyy=sum(Topel.Iyy(:))+sum(Botel.Iyy(:))+sum(Spars.Iyy(:));
Ixy=sum(Topel.Ixy(:))+sum(Botel.Ixy(:))+sum(Spars.Ixy(:));

format short e
disp('The Area Moment of Inertia are')
disp([Ixx Iyy Ixy])

%% Stress Analysis
% Constants
W=3200;         % m^2
E=20*10^6;      % Pa
g=9.8;          % kg*m/s^2
Wingspan=5.5;   % meters

% Set up
z=linspace(0,Wingspan,Wingspan*100+1);
Load=struct('L',zeros(6,length(z)), 'TotL',zeros(6,length(z)),...
            'D',zeros(6,length(z)), 'TotD',zeros(6,length(z)), ...
            'VL',zeros(6,length(z)), 'ML',zeros(6,length(z)),...
            'VD',zeros(6,length(z)), 'MD',zeros(6,length(z)),...
            'MM',zeros(6,length(z)), 'SigmaZ',zeros(6,length(z)));
        
% Critical Loads
% PHAA PLAA NHAA NLAA PosGust NegGust

% Speed
v=zeros(1,6);
% Load Factor
n=zeros(1,6);

% Calculate Lift Distribution
Load.L(1,:)=1;
LZ=zeros(1,length(z));

% Calculate Drag Distribution 
Load.D(1,:)=1;
DZ=zeros(1,length(z));

% Total Lift and Drag
for ii=1:length(z)-1
    Load.TotL(1,ii+1)=Load.TotL(1,ii)+...
                   0.5*(Load.L(1,ii)+Load.L(1,ii+1))*(z(ii+1)-z(ii));
    LZ(1,ii)=0.5*(Load.L(1,ii)+Load.L(1,ii+1))*(z(ii+1)-z(ii))*...
             (z(ii)+(z(ii+1)-z(ii))/2);
    Load.TotD(1,ii+1)=Load.TotD(1,ii)+...
                   0.5*(Load.D(1,ii)+Load.D(1,ii+1))*(z(ii+1)-z(ii));
    DZ(1,ii)=0.5*(Load.D(1,ii)+Load.D(1,ii+1))*(z(ii+1)-z(ii))*...
             (z(ii)+(z(ii+1)-z(ii))/2);
end
Lz=sum(LZ(1,:));
LZeq=Lz/Load.TotL(1,end);
Dz=sum(DZ(1,:));
DZeq=Dz/Load.TotD(1,end);

% Calculate Shear and Moment
Load.VL(1,1)= Load.TotL(1,end);
Load.ML(1,1)= Load.TotL(1,end)*LZeq;
Load.VD(1,1)= Load.TotD(1,end);
Load.MD(1,1)= Load.TotD(1,end)*DZeq;
ML=0;
MD=0;

for k=2:length(z)
    % Shear/Moment from Lift
    Load.VL(1,k)=Load.VL(1,1)-Load.TotL(1,k);
    ML=ML+0.5*(Load.VL(1,k-1)+Load.VL(1,k))*(z(k)-z(k-1));
    Load.ML(1,k)=Load.ML(1,1)-ML;
    
    % Shear/Moment from Drag
    Load.VD(1,k)=Load.VD(1,1)-Load.TotD(1,k);
    MD=MD+0.5*(Load.VD(1,k-1)+Load.VD(1,k))*(z(k)-z(k-1));
    Load.MD(1,k)=Load.MD(1,1)-MD;
end

% Moment
CM=-0.007;

% Check
figure()
hold on
plot(z,0,'k')
plot(z,Load.L(1,:),'r')
plot(z,Load.VL(1,:),'g')
plot(z,Load.ML(1,:),'b')
grid on
xlabel('Length (m)')
ylabel('Shear (N)')
%% Calculat Bending Stress and Deflection

% SigmaZ=(Mx*(Iyy*y-Ixy*x)+My*(Ixx*x-Ixy*y))/(Ixx*Iyy-Ixy^2);





