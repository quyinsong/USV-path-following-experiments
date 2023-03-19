function [Q_hat,z] = DiturbanceObserver( z, Hx, miu1, tau, Fx , B, ts)
%DITURBANCEOBSERVER �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
L = diag([5 5]);
Q_hat = L*(miu1-z);
z_dot = Fx+tau-B*Hx+Q_hat;
z = z+z_dot*ts;

end

