%% MAE 154B Wing Structure Analysis
close all;
clear all;
clc;
%% Finding Centroid
% All length are in unit of meters

skint=0.003;  % Skin thickness
spart=0.003; % Spar thickness
centroids=zeros(11,5);

% Creating nodes
nx=[0.02 0.08 0.30-spart/2 0.30+spart/2 0.914-spart/2 0.914+spart/2 1.35];
nyt=[0.03 0.06 0.1 0.1 0.075 0.075 0];
nyb=[-0.03 -.04 -.055 -.055 -.03 -.03 0];

% Airfoil
% Front
front=[nx(1)-2*skint nyb(1)-skint; nx(1)-2*skint nyt(1)+skint; nx(1) nyt(1)+skint; nx(1) nyb(1)-skint; nx(1)-2*skint nyb(1)-skint];
[centroids(1,1), centroids(1,2), centroids(1,3)]=CompCentroid(front(:,1),front(:,2));

% 1st Top
top1=[nx(1) nyt(1); nx(1) nyt(1)+skint; nx(2) nyt(2)+skint; nx(2) nyt(2); nx(1) nyt(1)];
[centroids(2,1), centroids(2,2), centroids(2,3)]=CompCentroid(top1(:,1),top1(:,2));
% 2nd Top
top2=[nx(2) nyt(2); nx(2) nyt(2)+skint; nx(3) nyt(3)+skint; nx(3) nyt(3); nx(2) nyt(2)];
[centroids(3,1), centroids(3,2), centroids(3,3)]=CompCentroid(top2(:,1),top2(:,2));
% 3rd Top
top3=[nx(4) nyt(4); nx(4) nyt(4)+skint; nx(5) nyt(5)+skint; nx(5) nyt(5); nx(4) nyt(4)];
[centroids(4,1), centroids(4,2), centroids(4,3)]=CompCentroid(top3(:,1),top3(:,2));
% 4th Top
top4=[nx(6) nyt(6); nx(6) nyt(6)+skint; nx(7) nyt(7)+skint; nx(7) nyt(7); nx(6) nyt(6)];
[centroids(5,1), centroids(5,2), centroids(5,3)]=CompCentroid(top4(:,1),top4(:,2));

% 1st Bottom
bot1=[nx(1) nyb(1); nx(1) nyb(1)-skint; nx(2) nyb(2)-skint; nx(2) nyb(2); nx(1) nyb(1)];
[centroids(6,1), centroids(6,2), centroids(6,3)]=CompCentroid(bot1(:,1),bot1(:,2));
% 2nd Bottom
bot2=[nx(2) nyb(2); nx(2) nyb(2)-skint; nx(3) nyb(3)-skint; nx(3) nyb(3); nx(2) nyb(2)];
[centroids(7,1), centroids(7,2), centroids(7,3)]=CompCentroid(bot2(:,1),bot2(:,2));
% 3rd Bottom
bot3=[nx(4) nyb(4); nx(4) nyb(4)-skint; nx(5) nyb(5)-skint; nx(5) nyb(5); nx(4) nyb(4)];
[centroids(8,1), centroids(8,2), centroids(8,3)]=CompCentroid(bot3(:,1),bot3(:,2));
% 4th Bottom
bot4=[nx(6) nyb(6); nx(6) nyb(6)-skint; nx(7) nyb(7)-skint; nx(7) nyb(7); nx(6) nyb(6)];
[centroids(9,1), centroids(9,2), centroids(9,3)]=CompCentroid(bot4(:,1),bot4(:,2));

% 1st Spar
spar1=[nx(3) nyb(3)-skint; nx(3) nyt(3)+skint; nx(4) nyt(4)+skint; nx(4) nyb(4)-skint; nx(3) nyb(3)-skint];
[centroids(10,1), centroids(10,2), centroids(10,3)]=CompCentroid(spar1(:,1),spar1(:,2));
% 2nd Spar
spar2=[nx(5) nyb(5)-skint; nx(5) nyt(5)+skint; nx(6) nyt(6)+skint; nx(6) nyb(6)-skint; nx(5) nyb(5)-skint];
[centroids(11,1), centroids(11,2), centroids(11,3)]=CompCentroid(spar2(:,1),spar2(:,2));

% Plot Wing Cross Section
figure()
hold on
plot(front(:,1),front(:,2))
plot(top1(:,1),top1(:,2))
plot(top2(:,1),top2(:,2))
plot(top3(:,1),top3(:,2))
plot(top4(:,1),top4(:,2))
plot(bot1(:,1),bot1(:,2))
plot(bot2(:,1),bot2(:,2))
plot(bot3(:,1),bot3(:,2))
plot(bot4(:,1),bot4(:,2))
plot(spar1(:,1),spar1(:,2))
plot(spar2(:,1),spar2(:,2))
xlim([-0.1,1.4])
ylim([-0.2,0.2])
xlabel('Length (m)')
ylabel('Length (m)')
grid on

% Calculating Centroid of the wing
TotalArea = sum(centroids(:,3));
for i=1:size(centroids,1)
    centroids(i,4)=centroids(i,1)*centroids(i,3);
    centroids(i,5)=centroids(i,2)*centroids(i,3);
    plot(centroids(i,1),centroids(i,2),'r*')
end
Cx=sum(centroids(:,4))/TotalArea;
Cy=sum(centroids(:,5))/TotalArea;

% Plot Centroid
plot(Cx,Cy,'b*')
hold off

%% Finding Moment of Inertia
Inertias=zeros(11,5);
[Inertias(1,1), Inertias(1,2), Inertias(1,3), Inertias(1,4), Inertias(1,5)]=CalcInertia(front(:,1),front(:,2),2*skint);
[Inertias(2,1), Inertias(2,2), Inertias(2,3), Inertias(2,4), Inertias(2,5)]=CalcInertia(top1(:,1),top1(:,2),skint);
[Inertias(3,1), Inertias(3,2), Inertias(3,3), Inertias(3,4), Inertias(3,5)]=CalcInertia(top2(:,1),top2(:,2),skint);
[Inertias(4,1), Inertias(4,2), Inertias(4,3), Inertias(4,4), Inertias(4,5)]=CalcInertia(top3(:,1),top3(:,2),skint);
[Inertias(5,1), Inertias(5,2), Inertias(5,3), Inertias(5,4), Inertias(5,5)]=CalcInertia(top4(:,1),top4(:,2),skint);
[Inertias(6,1), Inertias(6,2), Inertias(6,3), Inertias(6,4), Inertias(6,5)]=CalcInertia(bot1(:,1),bot1(:,2),skint);
[Inertias(7,1), Inertias(7,2), Inertias(7,3), Inertias(7,4), Inertias(7,5)]=CalcInertia(bot2(:,1),bot2(:,2),skint);
[Inertias(8,1), Inertias(8,2), Inertias(8,3), Inertias(8,4), Inertias(8,5)]=CalcInertia(bot3(:,1),bot3(:,2),skint);
[Inertias(9,1), Inertias(9,2), Inertias(9,3), Inertias(9,4), Inertias(9,5)]=CalcInertia(bot4(:,1),bot4(:,2),skint);
[Inertias(10,1), Inertias(10,2), Inertias(10,3), Inertias(10,4), Inertias(10,5)]=CalcInertia(spar1(:,1),spar1(:,2),spart);
[Inertias(11,1), Inertias(11,2), Inertias(11,3), Inertias(11,4), Inertias(11,5)]=CalcInertia(spar2(:,1),spar2(:,2),spart);

% Apply Parallel Axis 
for k=1:size(Inertias,1)
     Inertias(k,3)=Inertias(k,3)+centroids(k,3)*(centroids(k,2)-Cy)^2;
     Inertias(k,4)=Inertias(k,4)+centroids(k,3)*(centroids(k,1)-Cx)^2;
     Inertias(k,5)=Inertias(k,5)+centroids(k,3)*(Cx-centroids(k,1))*(Cy-centroids(k,2));
end

TotalIxx=sum(Inertias(:,3));
TotalIyy=sum(Inertias(:,4));
TotalIxy=sum(Inertias(:,5));

%% Output
format short e
disp('The centroid is at')
disp([Cx,Cy])
disp('The Area Moment of Inertia are')
disp([TotalIxx TotalIyy TotalIxy])

%% Stress Analysis







