function [x, y, A]=CompCentroid(Xarray, Yarray)
x=0.5*((Xarray(1)+Xarray(2))/2+(Xarray(3)+Xarray(4))/2);
y=0.5*((Yarray(1)+Yarray(2))/2+(Yarray(3)+Yarray(4))/2);

a=[Xarray(1)-Xarray(2); Yarray(1)-Yarray(2); 0];
b=[Xarray(3)-Xarray(2); Yarray(3)-Yarray(2); 0];
A=sum(abs(cross(a,b)));
end
