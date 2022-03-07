function xnext = USV03( x,tao,d )
%USV03 
%   Author: Quyinsong
%   Date: 5th March 2022 
%   Reference: 1. 固定时间预测器下的欠驱动无人艇路径跟踪控制
% 2. Planar trajectory planning and tracking control design for underactuated AUVs
ts = 0.01;
%% USV03 parameters
m11 = 215; m22 = 265; m33 = 80;
Xu = 70; Xuu = 100; Yv = 100; Yvv = 200; Nr = 50; Nrr = 100;
%% USV03 state
u = x(1); v = x(2); r = x(3); psi = x(6);
Fu = tao(1); Tr = tao(3);
V = [u v r]';
Rbn = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
%% time derevatives
Xdot = Rbn*V;
Vdot = [m22*v*r/m11-Xu*u/m11-Xuu*abs(u)*u+Fu/m11+d(1);
        -m11*u*r/m22-Yv*v/m22-Yvv*abs(v)*v/m22+d(2);
        (m11-m22)*u*v/m33-Nr*r/m33-Nrr*abs(r)*r+Tr/m33+d(3)];
xdot = [Vdot;Xdot];
%% state update
xnext = euler2(xdot,x,ts);

end

