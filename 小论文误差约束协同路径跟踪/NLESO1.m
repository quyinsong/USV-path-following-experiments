function [xe1_hat, ye1_hat, u_hat,v_hat ] = NLESO1( xe,ye,kc,up,psi,psip,ts )
%LESO 此处显示有关此函数的摘要
%   xe: 纵向跟踪误差
%   ye: 横向跟踪误差
%   kc：曲率
%   up：路径虚拟点速度
%   ts：采样时间
persistent xe_hat ye_hat g1_hat g2_hat kx1 kx2 ky1 ky2 rx ry;
if isempty(xe_hat)
    xe_hat = xe;
    ye_hat = ye;
    g1_hat = 0;
    g2_hat = 0;
    kx1 = 1;
    kx2 = 1;
    ky1 = 1;
    ky2 = 1;
    rx = 5; ry = 5;
end



xe_hat_dot = g1_hat+kc*up*ye-up-kx1*fal(rx*(xe_hat-xe), 0.7, 1);
g1_hat_dot = -rx*kx2*fal(rx*(xe_hat-xe),0.4,1);
xe_hat = xe_hat_dot*ts+xe_hat;
g1_hat = g1_hat_dot*ts+g1_hat;

ye_hat_dot = g2_hat-kc*up*xe-ky1*fal(ry*(ye_hat-ye),0.7,1);
g2_hat_dot = -ry*ky2*fal(ry*(ye_hat-ye),0.4,1);
ye_hat = ye_hat_dot*ts+ye_hat;
g2_hat = g2_hat_dot*ts+g2_hat;

Rp = [cos(psi-psip) -sin(psi-psip);
      sin(psi-psip) cos(psi-psip)];
  
eta2_hat = Rp\[g1_hat;g2_hat];

u_hat = g1_hat*cos(psi-psip)+g2_hat*sin(psi-psip);
v_hat = g2_hat*cos(psi-psip)-g1_hat*sin(psi-psip);

xe1_hat = xe_hat;
ye1_hat = ye_hat;
end


