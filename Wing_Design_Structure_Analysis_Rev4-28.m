%% MAE 154B Wing Design Structure Analysis
close all;
clear all;
clc;

%% Define Nodes
% Import Shape of NACA2412
TopNodes = importdata('NACA2412_Top.txt');
BottomNodes = importdata('NACA2412_Bottom.txt');
% Extract Node Points
NX=TopNodes(:,1);
NYtop=TopNodes(:,2);
NYbot=BottomNodes(:,2);
% Scale points to desire size
scalefactor=1.346;
nx=scalefactor*NX;
nytop=scalefactor*NYtop;
nybot=scalefactor*NYbot;
nx=nx(1:15);
% Fitting the shape to a polynomial
% [Ptop, Stop]=polyfit(Nx,Nytop,12);
% [Pbot, Sbot]=polyfit(Nx,Nybot,12);
% Set nodes every 2 mm for first 80% of cord length
% nx=0:0.002:0.85*1.346;
% [nytop, topdelta]=polyval(Ptop,nx,Stop);
% [nybot, botdelta]=polyval(Pbot,nx,Sbot);
% Test fit
% figure()
% hold on
% plot(Nx,Nytop,'r')
% plot(nx,nytop)
% plot(Nx,Nybot,'r')
% plot(nx,nybot)
%% Create Elements
% All Nodes for airfoil are in the inner side of the elements
skint=0.002;  % Skin thickness
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
%    plot(Topel.cx(i), Topel.cy(i),'r*')
   
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
%    plot(Botel.cx(i), Botel.cy(i),'r*')
   
   XiAiairf=XiAiairf+Topel.cx(i)*Topel.Area(i)+Botel.cx(i)*Botel.Area(i);
   YiAiairf=YiAiairf+Topel.cy(i)*Topel.Area(i)+Botel.cy(i)*Botel.Area(i);
end
% Centroid of Airfoil
AirfoilArea=sum(Topel.Area(:,1))+sum(Botel.Area(:,1));
Cxaf=XiAiairf/AirfoilArea;
Cyaf=YiAiairf/AirfoilArea;

% Create Spars/Spar-Caps Elements
spart=0.004; % Spar thickness
sparcapt=0.002; % Spar cap size

% Spar Location
SparIndex=[find(nx>0.3*1.346-0.001 & nx<0.3*1.346+0.001) length(nx)];
Spars=struct('posX',nx(SparIndex),...
             'posY',nybot(SparIndex),...
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

%Spar Caps
%Assume Spar caps are angle with equal length legs
SparCaps=struct('posX',transpose(Spars.xcoord(:,1:4)),...
                'posY',transpose(Spars.ycoord(:,1:4)),...
                'Area',zeros(6,1),...
                'cx',zeros(6,1), ...
                'cy',zeros(6,1),...
                'Ixx',zeros(6,1),...
                'Iyy',zeros(6,1),...
                'Ixy',zeros(6,1),...
                'xcoord',zeros(6,7),...
                'ycoord',zeros(6,7));  
b=2*sparcapt;
XiAisparcap=0;
YiAisparcap=0;
for i=1:6
    % Area of spar caps
    SparCaps.Area(i)=(2*b-sparcapt)*sparcapt;
    % centroids of spar caps
    if i==3||i==4
        SparCaps.cx(i)=SparCaps.posX(i)+(1/SparCaps.Area(i))*...
                       (sparcapt/2*(b^2+b*sparcapt-sparcapt^2));
    else
        SparCaps.cx(i)=SparCaps.posX(i)-(1/SparCaps.Area(i))*...
                        (sparcapt/2*(b^2+b*sparcapt-sparcapt^2));
    end
    
    if i==1||i==4||i==5
        SparCaps.cy(i)=SparCaps.posY(i)+(1/SparCaps.Area(i))*...
                       (sparcapt/2*(b^2+b*sparcapt-sparcapt^2));
    else
        SparCaps.cy(i)=SparCaps.posY(i)-(1/SparCaps.Area(i))*...
                        (sparcapt/2*(b^2+b*sparcapt-sparcapt^2));
    end
    % Moment of inertia
    SparCaps.Ixx(i)=(sparcapt/3)*(b*sparcapt^2+b^3-sparcapt^3)-...
                     SparCaps.Area(i)*(SparCaps.cy(i)-SparCaps.posY(i))^2;
    SparCaps.Iyy(i)=(sparcapt/3)*(b*sparcapt^2+b^3-sparcapt^3)-...
                     SparCaps.Area(i)*(SparCaps.cx(i)-SparCaps.posX(i))^2;
    SparCaps.Ixy(i)=(sparcapt^2/4)*(2*b^2-sparcapt^2)-SparCaps.Area(i)*...
                    (SparCaps.cy(i)-SparCaps.posY(i))*...
                    (SparCaps.cx(i)-SparCaps.posX(i));
                
    XiAisparcap=XiAisparcap+SparCaps.cx(i)*SparCaps.Area(i);
    YiAisparcap=YiAisparcap+SparCaps.cy(i)*SparCaps.Area(i);
    % Draw Spar Caps
    if i==3||i==4
        SparCaps.xcoord(i,:)=[SparCaps.posX(i) SparCaps.posX(i)+b...
                          SparCaps.posX(i)+b SparCaps.posX(i)+sparcapt...
                          SparCaps.posX(i)+sparcapt SparCaps.posX(i)...
                          SparCaps.posX(i)];
    else
        SparCaps.xcoord(i,:)=[SparCaps.posX(i) SparCaps.posX(i)-b...
                          SparCaps.posX(i)-b SparCaps.posX(i)-sparcapt...
                          SparCaps.posX(i)-sparcapt SparCaps.posX(i)...
                          SparCaps.posX(i)];
    end

    if i==1||i==4||i==5
        SparCaps.ycoord(i,:)=[SparCaps.posY(i), SparCaps.posY(i),...
                          SparCaps.posY(i)+sparcapt,...
                          SparCaps.posY(i)+sparcapt,...
                          SparCaps.posY(i)+b, SparCaps.posY(i)+b,...
                          SparCaps.posY(i)];
    else
        SparCaps.ycoord(i,:)=[SparCaps.posY(i) SparCaps.posY(i)...
                          SparCaps.posY(i)-sparcapt...
                          SparCaps.posY(i)-sparcapt...
                          SparCaps.posY(i)-b SparCaps.posY(i)-b...
                          SparCaps.posY(i)];
    end
    
    plot(SparCaps.xcoord(i,:),SparCaps.ycoord(i,:))
    plot(SparCaps.cx(i),SparCaps.cy(i),'r*')
end
% Centroid of SparCaps
SparCapArea=sum(SparCaps.Area(:,1));
Cxsparcap=XiAisparcap/SparCapArea;
Cysparcap=YiAisparcap/SparCapArea;

% Calculate Total Centroid
Cx=(XiAiairf+XiAiSpar+XiAisparcap)/(AirfoilArea+SparArea+SparCapArea);
Cy=(YiAiairf+YiAiSpar+YiAisparcap)/(AirfoilArea+SparArea+SparCapArea);

plot(Cx,Cy,'r*')
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

Ixx=sum(Topel.Ixx(:))+sum(Botel.Ixx(:))+sum(Spars.Ixx(:))+...
        sum(SparCaps.Ixx(:));
Iyy=sum(Topel.Iyy(:))+sum(Botel.Iyy(:))+sum(Spars.Iyy(:))+...
        sum(SparCaps.Iyy(:));
Ixy=sum(Topel.Ixy(:))+sum(Botel.Ixy(:))+sum(Spars.Ixy(:))+...
        sum(SparCaps.Ixy(:));

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
Load=struct('Wy',zeros(6,length(z)), 'TotWy',zeros(6,length(z)),...
            'Wx',zeros(6,length(z)), 'TotWx',zeros(6,length(z)), ...
            'VWy',zeros(6,length(z)), 'MWy',zeros(6,length(z)),...
            'VWx',zeros(6,length(z)), 'MWx',zeros(6,length(z)),...
            'MM',zeros(6,length(z)), 'SigmaZ',zeros(6,length(z)));
        
% Critical Loads
% PHAA PLAA NHAA NLAA PosGust NegGust

% Speed
v=zeros(1,6);
% Load Factor
n=zeros(1,6);

% Calculate Lift Distribution
Load.Wy(1,:)=1;
WyZ=zeros(1,length(z));

% Calculate Drag Distribution 
Load.Wx(1,:)=1;
WxZ=zeros(1,length(z));

% Total Lift and Drag
for ii=1:length(z)-1
    Load.TotWy(1,ii+1)=Load.TotWy(1,ii)+...
                   0.5*(Load.Wy(1,ii)+Load.Wy(1,ii+1))*(z(ii+1)-z(ii));
    WyZ(1,ii)=0.5*(Load.Wy(1,ii)+Load.Wy(1,ii+1))*(z(ii+1)-z(ii))*...
             (z(ii)+(z(ii+1)-z(ii))/2);
    Load.TotWx(1,ii+1)=Load.TotWx(1,ii)+...
                   0.5*(Load.Wx(1,ii)+Load.Wx(1,ii+1))*(z(ii+1)-z(ii));
    WxZ(1,ii)=0.5*(Load.Wx(1,ii)+Load.Wx(1,ii+1))*(z(ii+1)-z(ii))*...
             (z(ii)+(z(ii+1)-z(ii))/2);
end
Wyz=sum(WyZ(1,:));
WyZeq=Wyz/Load.TotWy(1,end);
Wxz=sum(WxZ(1,:));
WxZeq=Wxz/Load.TotWx(1,end);

% Calculate Shear and Moment
Load.VWy(1,1)= Load.TotWy(1,end);
Load.MWy(1,1)= Load.TotWy(1,end)*WyZeq;
Load.VWx(1,1)= Load.TotWx(1,end);
Load.MWx(1,1)= Load.TotWx(1,end)*WxZeq;
MWy=0;
MWx=0;

for k=2:length(z)
    % Shear/Moment from Lift
    Load.VWy(1,k)=Load.VWy(1,1)-Load.TotWy(1,k);
    MWy=MWy+0.5*(Load.VWy(1,k-1)+Load.VWy(1,k))*(z(k)-z(k-1));
    Load.MWy(1,k)=Load.MWy(1,1)-MWy;
    
    % Shear/Moment from Drag
    Load.VWx(1,k)=Load.VWx(1,1)-Load.TotWx(1,k);
    MWx=MWx+0.5*(Load.VWx(1,k-1)+Load.VWx(1,k))*(z(k)-z(k-1));
    Load.MWx(1,k)=Load.MWx(1,1)-MWx;
end

% Moment
CM=-0.007;

% Check
figure()
hold on
plot(z,zeros(1,length(z)),'k', 'Linewidth',3)
plot(z,Load.Wy(1,:),'r')
hAx=plotyy(z,Load.VWy(1,:),z,Load.MWy(1,:));
legend('Beam','Applied Force','Shear','Moment')
xlabel('Length (m)')
ylabel(hAx(1), 'Shear (N)')
ylabel(hAx(2), 'Moment (N-m)')
grid on
hold off
%% Calculat Bending Stress and Deflection
% SigmaZ=(Mx*(Iyy*y-Ixy*x)+My*(Ixx*x-Ixy*y))/(Ixx*Iyy-Ixy^2);





