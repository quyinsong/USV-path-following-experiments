function [xe1_hat, ye1_hat, u_hat,v_hat ] = LESO1( xe,ye,kc,up,psi,psip,ts )
%LESO �˴���ʾ�йش˺�����ժҪ
%   xe: ����������
%   ye: ����������
%   kc������
%   up��·��������ٶ�
%   ts������ʱ��
persistent xe_hat ye_hat g1_hat g2_hat kx1 kx2 ky1 ky2;
if isempty(xe_hat)
    xe_hat = xe;
    ye_hat = ye;
    g1_hat = 0;
    g2_hat = 0;
    kx1 = 1;
    kx2 = 1;
    ky1 = 1;
    ky2 = 1;
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

u_hat = g1_hat*cos(psi-psip)+g2_hat*sin(psi-psip);
v_hat = g2_hat*cos(psi-psip)-g1_hat*sin(psi-psip);

xe1_hat = xe_hat;
ye1_hat = ye_hat;
end

