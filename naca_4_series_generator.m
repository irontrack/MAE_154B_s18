c = 5;
m = 0.02;
p = 0.4;
t = 0.12;
pc = p*c;
x = linspace(0,c,500);

yt = 5*t*(0.2969.*sqrt(x/c) - 0.1260.*(x/c) - 0.3516.*(x/c).^2 + 0.2843.*(x/c).^3 - 0.1015.*(x/c).^4);
yc = zeros(1,500);
for i = 1:500
    if x(i) <= pc
        yc(i) = (m/(p^2))*(2*p*(x(i)/c) - (x(i)/c)^2);
    else
        yc(i) = (m/(1 - p)^2)*((1 - 2*p) + 2*p*(x(i)/c) - (x(i)/c)^2);
    end
end

dycdx = zeros(1,500);
for i = 1:500
    if x <= pc
        dycdx(i) = (2*m/p^2)*(p - (x(i)/c));
    else
        dycdx(i) = (2*m/(1 - p)^2)*(p - (x(i)/c));
    end
end
theta = atan(dycdx);

XU = x - yt.*sin(theta);
XL = x + yt.*sin(theta);
YU = yc + yt.*cos(theta);
YL = yc - yt.*cos(theta);

plot(XU,YU)
hold on
plot(XL,YL)
grid
axis([0,c,-t*c,t*c])