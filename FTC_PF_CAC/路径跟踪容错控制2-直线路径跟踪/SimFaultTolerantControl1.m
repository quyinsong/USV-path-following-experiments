% test fault-tolerant control
% date: 7/2/2022
% Author: Quyinsong
clc
clear all
close all

ts = 0.02;
tfinal = 80;
Ns = tfinal/ts;
tf = 30;
x = [0 0 0 8 0 pi/3]'; w = 0;
Tmax = 113; Tmin = -55;
% fault detection
D = [0;0];
L = diag([20 20]);
D_up = 0.01;
Tc1 = [0 0]';
tauw = 0;
% adaptive RBFNN fault observer 
Q = [0 0]'; Hx_hat = [0 0]'; Hx = [0 0]'; Q_hat = [0 0]'; z = 0;
% Auxliary system
A = diag([50 50]);
lamda = [0 0]';

m11 = 50.05; m33 = 17.21; Fx = [0 0]';
B = [1/m11 1/m11; 0.26/m33 -0.26/m33];
% LOS
pc = 1; 
tau = [0 0]';
k1 = 3; delta = 2;  T1 = 0.1; T2 = 0.1; T3 = 0.2;
k2 = 0.8;  x1 = 0; x2 = 0; ud = 0.6;
x3 = [0;0]; x3_d = [0;0];
% control law
K = diag([1 1]);

for k = 1:Ns
    t = (k-1)*ts;
    % diturbances
    switch tauw
        case 1
            tau_w = [1.2*sin(t)*cos(0.5*t) 0.5*sin(t)*cos(0.5*t) 1.1*sin(t)*cos(0.5*t)]';
        case 0
            tau_w = [0 0 0]';
    end
    
    % state
    xn = x(4); yn = x(5); psi = x(6); u = x(1); v = x(2); r = x(3); U = sqrt(u^2+v^2); beta = atan2(v,u);
    % LOS guidance
    [xd, yd, Xsf_d, psif, kc, up , xe, ye, w] = LOS( pc, k1, delta, xn, yn, psi, U, beta, w, ts );
    % backstepping control law (virtual control law)
    [ x1, x1_dot] = TD1( beta-Xsf_d,x1,T1,ts );
    Xsf = psi+beta-psif;
    Xsf_tilde = ssa(Xsf-Xsf_d);
    alpha_r = -x1_dot+kc*up-k2*Xsf_tilde;
    [ x2, x2_dot] = TD1( alpha_r,x2,T2,ts );
    e_r = r - x2;
    e_u = u-ud;
    
    % fault detection and diagnosis
    miu1 = [u,r]';
%     Fx = [1.6855*v*r-3.1483*u
%           -1.9936*u*v-2.0081*r];
    
    [D,D_tilde] = FaultDetection( D, Fx, miu1, B*Tc1, ts, L);
    % decision mechanism
    if norm(D_tilde) >= D_up
        flag = 1;
    else
        flag = 0;
    end
    % adaptive RBFNN fault observer       
    Ed_dot = [0,x2_dot]';
    E = [e_u,e_r]'-lamda;
    switch flag
        case 1
            [Q, Hx_hat, D_hat ] = FaultEstimator1( miu1,B*Tc1,B,ts, Fx );
        case 0
            Hx_hat = [0 0]';
    end
    
    if flag == 1
        tau = -Fx-K*E+Ed_dot-A*lamda-0.05*sign(E)+B*Hx;
    else
        tau = -Fx-K*E+Ed_dot-A*lamda-0.05*sign(E);
    end
    
    Tc = B\tau;
    Tc0 = Tc;
    % 执行器饱和
    for i = 1:2
        if Tc(i)>=Tmax
            Tc(i) = Tmax;
        elseif Tc(i)<=Tmin
            Tc(i)=Tmin;
        end
    end
    Tc1 = Tc;
    % Auxliary system
    lamda_dot = -A*lamda-B*(Tc0-Tc1);
    lamda = euler2(lamda_dot,lamda,ts);
    % 惯性环节
    for i=1:2
        [x3(i),x3_d(i)] = TD1(Tc1(i),x3(i),T3,ts);
        Tc1(i) = x3(i);
    end
    % USV
    [x, Tf, Hx, Fx] = USV( x, Tc1, tau_w, ts, t, tf );
    % estimate value
    xout(k,:) = [t x' Tc1' Tf' xd yd xe ye D_tilde' norm(D_tilde) flag Q' Hx_hat' Hx'];
end
t = xout(:,1);
u = xout(:,2);
v = xout(:,3);
r = xout(:,4);
xn = xout(:,5);
yn = xout(:,6);
psi = xout(:,7);
Tc1 = xout(:,8:9);
Tf = xout(:,10:11);
xd = xout(:,12);
yd = xout(:,13);
xe = xout(:,14);
ye = xout(:,15);
D_tilde = xout(:,16:17);
norm_Dtilde = xout(:,18);
flag = xout(:,19);
miu1_hat = xout(:,20:21);
Hx_hat = xout(:,22:23);
Hx = xout(:,24:25);

figure(1);hold on
xrange=[-10 3 30]; yrange = [-10 3 30];
for k=1:1:Ns
    pos =[xn(k) yn(k)]';
    if k==1
        modelplot(pos,psi(k),xrange,yrange);
    end
    if rem(k,500)==0
        modelplot(pos,psi(k),xrange,yrange);
    end   
end
Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);
plot(yn,xn,'r-','linewid',1);
plot(yd,xd,'b-','linewid',1);
set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
axis([Xmin Xmax,Ymin Ymax]);
grid
hold off;

figure(2)
plot(t,u,'r-','linewid',0.5);
figure(3)
plot(t,xe,'r-',t,ye,'b-','linewid',1);
legend('xe','ye');
figure(4)
plot(t,D_tilde(:,1),'r-',t,D_tilde(:,2),'b-',t,norm_Dtilde,'k-')
figure(5)
plot(t,Tc1(:,1),'r-',t,Tc1(:,2),'b-',t,Tf(:,1),'r--',t,Tf(:,2),'b--')
figure(6)
plot(t,miu1_hat(:,1),'r-',t,miu1_hat(:,2),'b-'); hold on
plot(t,u,'r--',t,r,'b--'); hold off
legend('u_hat','r_hat','u','v');
figure(7)
plot(t,flag);
figure(8)
plot(t,-Hx_hat(:,1),'r-',t,-Hx_hat(:,2),'b-'); hold on
plot(t,-Hx(:,1),'r--',t,-Hx(:,2),'b--'); hold off
legend('h1hat','h2hat','h1','h2');
