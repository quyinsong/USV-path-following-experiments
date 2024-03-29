function [xout, xdot, T, fu, fr] = ASV3( x0, Tc, tau_w, Vc, ts)
%USV USV model 
% ControlInput: T1 and T2
% DiturbanceInput: tau_w
% SampleTimeInput: ts
% SimulationTimeInput: t
% FaultTime: tf
% State: x=[xn yn psi u v r]';
% Output: x
% NumericalSimulationMethod: 1-nd Euler 
% Reference: Modeling and Experimental Testing of an Unmanned
%            Surface Vehicle with Rudderless Double Thrusters
% Date: 6/30/2022/
% Author: Quyinsong
% 
persistent x
if isempty(x)
   x = x0; 
end
%% State
xn = x(4); yn = x(5); psi = x(6); 
u = x(1); v = x(2); r = x(3);
if abs(r)>=0.5
    r = sign(r)*0.5;
end
x = [u,v,r,xn,yn,psi]';
miu = [x(1) x(2) x(3)]'; 
%% USV model parameters known perfectly
m11 = 50.5; m22 = 84.36; m33 = 17.21;
d11 = 151.57; d22 = 132.5; d33 = 34.56;
C = [0     0 -m22*v
     0     0  m11*u
     m22*v -m11*u 0];
D = diag([d11 d22 d33]);
M = diag([m11 m22 m33]);

Umax = 1.5;

%% thruster
dp = 0.26;
% B = [1 1;0 0; dp -dp];
B = [1/m11 1/m11;0 0; dp/m33 -dp/m33];

Tmax = Umax*d11/2; Tmin = -55;
Tdot_max = 100; Tdot_min = -100;
persistent Tf Tt
if isempty(Tf)
    Tf = Tc;
    Tt = 0.5; % 执行器积分常数
end

Tf_dot = -(Tf-Tc)/Tt;

Tf = Tc;

% 执行器速率限幅
for i = 1:2
    if Tf_dot(i)>=Tdot_max
        Tf_dot(i) = Tdot_max;
    elseif Tf_dot(i)<=Tdot_min
        Tf_dot(i) = Tdot_min;
    end
end

Tf = euler2(Tf_dot,Tf,ts);

% 执行器输出限幅
for i = 1:2
    if Tf(i)>=Tmax
        Tf(i) = Tmax;
    elseif Tf(i)<=Tmin
        Tf(i)=Tmin;
    end
end

% tau = B*Tf;

%% state update
R = [cos(psi) -sin(psi) 0
     sin(psi) cos(psi)  0
     0        0         1];
F = M\(-D*miu-C*miu+tau_w);
fu = F(1); fv = F(2); fr = F(3);

% miu_dot = M\(-D*miu-C*miu+tau+tau_w);
miu_dot = F+B*Tf;
eta_dot = R*miu+[Vc;0];

xdot = [miu_dot;eta_dot];
x = euler2(xdot,x,ts);
xout = x;
% 饱和输出
T = Tf;
end

