function [ x, vd_d ] = TD1( vd,x,T,ts )
%DSC 1-st TD

x_d = -(x-vd)/T;
x = euler2(x_d,x,ts);
vd_d = x_d;
end

