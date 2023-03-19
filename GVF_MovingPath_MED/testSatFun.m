% 测试非连续饱和函数的连续近似
clc
clear all
close all

x = -10:0.1:10;
x_up = 2; x_dn = -1;
for i = 1:length(x)
    if x(i)<x_dn
        y(i) = x_dn;
    elseif x(i)<=x_up
        y(i) = x(i);
    else
        y(i) = x_up;
    end
    if x(i)<0
        y1(i) = x_dn*tanh(x(i)/x_dn);
    else
        y1(i) = x_up*tanh(x(i)/x_up);
    end
end

figure(1)
plot(x,y,'r-',x,y1,'b-','linewid',2);