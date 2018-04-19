function [length, angle, Ixx, Iyy, Ixy]=CalcInertia(Xarray, Yarray, t)
length=sqrt((Xarray(3)-Xarray(2))^2+(Yarray(3)-Yarray(2))^2);
angle=atan2((Yarray(3)-Yarray(2)),(Xarray(3)-Xarray(2)));
    if angle==0
        Ixx=t*length^3/12;
        Iyy=t^3*length/12;
        Ixy=0;
    else
        Ixx=length^3*t*(sin(angle))^2/12;
        Iyy=length^3*t*(cos(angle))^2/12;
        Ixy=length^3*t*sin(2*angle)/24;
    end
end