% test LOS backstepping
% date: 6/30/2022
% Author: Quyinsong
clc
clear all
close all

ts = 0.02;
tfinal = 100;
Ns = tfinal/ts;
tf = 50;
x = [0 0 0 5 0 0]'; w = 0;
tau_w = [0 0 0]';
B = [1 1; 0.26 -0.26];
pc = 1; k1 = 3; delta = 3;  T1 = 0.1; T2 = 0.2; T3 = 0.5;
k2 = 0.5; k3 = 1; k4 = 0.8; x1 = 0; x2 = 0; 
x3 = [0;0]; x3_d = [0;0];
for k = 1:Ns
    t = (k-1)*ts;
    tau_w = [5*sin(t) 2*sin(t) 2*sin(t)]';
    xn = x(4); yn = x(5); psi = x(6); u = x(1); v = x(2); r = x(3); U = sqrt(u^2+v^2); beta = atan2(v,u);
    [xd, yd, alpha_psie, psif, psif_dot, xe, ye, w] = LOS( pc, k1, delta, xn, yn, psi, U, beta, w, ts );
    [ x1, x1_dot] = TD1( alpha_psie,x1,T1,ts );
    e_psi = ssa(psi-psif-x1);
    alpha_r = psif_dot+x1_dot-k2*e_psi;
    [ x2, x2_dot] = TD1( alpha_r,x2,T2,ts );
    e_r = r - x2;
    tau = [50.05*(3.1483*u-1.6855*v*r-k4*(u-0.8))
           17.21*(1.9936*u*v+2.0081*r+x2_dot-k3*e_r)];

    T = B\tau;
    for i=1:2
        [x3(i),x3_d(i)] = TD1(T(i),x3(i),T3,ts);
        T(i) = x3(i);
    end
    [x, Tf] = USV( x, T, tau_w, ts, t, tf );
    
    xout(k,:) = [t x' T' Tf' xd yd xe ye];
end
t = xout(:,1);
u = xout(:,2);
v = xout(:,3);
r = xout(:,4);
xn = xout(:,5);
yn = xout(:,6);
psi = xout(:,7);
T = xout(:,8:9);
Tf = xout(:,10:11);
xd = xout(:,12);
yd = xout(:,13);
xe = xout(:,14);
ye = xout(:,15);

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
plot(t,u,'linewid',2);
figure(3)
plot(t,xe,'r-',t,ye,'b-','linewid',1);
legend('xe','ye');
figure(4)
plot(t,T(:,1),'r-',t,T(:,2),'b-','linewid',1); hold on
plot(t,Tf(:,1),'r--',t,Tf(:,2),'b--','linewid',1); hold off
legend('T1','T2','Tf1','Tf2');