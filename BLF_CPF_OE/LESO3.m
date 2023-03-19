function [xe1_hat, ye1_hat, u_hat,v_hat ] = LESO3( xe,ye,kc,up,psi,psip,ts )
%LESO 此处显示有关此函数的摘要
%   xe: 纵向跟踪误差
%   ye: 横向跟踪误差
%   kc：曲率
%   up：路径虚拟点速度
%   ts：采样时间
persistent xe_hat ye_hat g1_hat g2_hat kx1 kx2 ky1 ky2;
if isempty(xe_hat)
    xe_hat = xe;
    ye_hat = ye;
    g1_hat = 0;
    g2_hat = 0;
    kx1 = 5;
    kx2 = 8;
    ky1 = 5;
    ky2 = 8;
end

xe_hat_dot = g1_hat+kc*up*ye-up-kx1*(xe_hat-xe);
g1_hat_dot = -kx2*(xe_hat-xe);
xe_hat = xe_hat_dot*ts+xe_hat;
g1_hat = g1_hat_dot*ts+g1_hat;

ye_hat_dot = g2_hat-kc*up*xe-ky1*(ye_hat-ye);
g2_hat_dot = -ky2*(ye_hat-ye);
ye_hat = ye_hat_dot*ts+ye_hat;
g2_hat = g2_hat_dot*ts+g2_hat;

Rp = [cos(psi-psip) -sin(psi-psip);
      sin(psi-psip) cos(psi-psip)];
  
eta2_hat = Rp\[g1_hat;g2_hat];

u_hat = eta2_hat(1);
v_hat = eta2_hat(2);
xe1_hat = xe_hat;
ye1_hat = ye_hat;
end

