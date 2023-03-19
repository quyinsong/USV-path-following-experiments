function [Q_hat,z] = DiturbanceObserver( z, Hx, miu1, tau, Fx , B, ts)
%DITURBANCEOBSERVER 此处显示有关此函数的摘要
%   此处显示详细说明
L = diag([5 5]);
Q_hat = L*(miu1-z);
z_dot = Fx+tau-B*Hx+Q_hat;
z = z+z_dot*ts;

end

