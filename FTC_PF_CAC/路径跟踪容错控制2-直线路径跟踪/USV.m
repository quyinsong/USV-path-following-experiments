function [x, Tf, Hx, Fx] = USV( x, T, tau_w, ts, t, tf )
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
%% State
xn = x(4); yn = x(5); psi = x(6); eta = [xn yn psi]';
u = x(1); v = x(2); r = x(3); miu = [u v r]';
%% USV model parameters known perfectly
m11 = 50.5; m22 = 84.36; m33 = 17.21;
d11 = 151.57; d22 = 132.5; d33 = 34.56;
C = [0     0 -m22*v
     0     0  m11*u
     m22*v -m11*u 0];
D = diag([d11 d22 d33]);
M = diag([m11 m22 m33]);
Fx3 = M\(-C*miu-D*miu);
Fx = [Fx3(1),Fx3(3)]';
Umax = 1.5;
%% USV model unkown terms
DeltaFx = [0.02*d11*abs(u)*u; 0.05*d22*abs(v)*v; 0.02*d33*abs(r)*r];
%% thruster
dp = 0.26;
B = [1 1;0 0; dp -dp];
a = [0 0]';
beta = [0 0]';
zeta = [0 0]';
Tf = T;
rho = zeros(2,1);
Tmax = Umax*d11/2; Tmin = -55;
Hx = zeros(2,1);
for i = 1:2
    if t>=15 && t<35
        a = [1 1]'; beta = [2 2]'; zeta = [0 0]';
    end
    if t>=35 && t<55
        a = [1 1]'; beta = [2 99]'; zeta = [0 0]';
    end
    if t>=55 && t<75
        a = [1 1]'; beta = [2 99]'; zeta = [0.05*sin(0.8*t) 0.05*cos(0.8*t)]';
    end
    if t>=75
        a = [0 0]';
    end
    rho(i)=a(i)/(1+beta(i)*exp(-zeta(i)));
    Tf(i) = (1-rho(i))*T(i);
    Hx(i) = rho(i)*T(i);
    if Tf(i)>=Tmax
        Tf(i) = Tmax;
    elseif Tf(i)<=Tmin
        Tf(i)=Tmin;
    end
end

DeltaB = 0.01*B;
tau = (B+DeltaB)*Tf;

%% state update
R = [cos(psi) -sin(psi) 0
     sin(psi) cos(psi)  0
     0        0         1];
miu_dot = M\(-D*miu-C*miu+tau+tau_w);
eta_dot = R*miu;

xdot = [miu_dot;eta_dot];
x = euler2(xdot,x,ts);


end

