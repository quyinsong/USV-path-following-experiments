% test USV 
% date: 6/30/2022
% Author: Quyinsong
clc
clear all
close all

ts = 0.02;
tfinal = 80;
Ns = tfinal/ts;
tf = 40;
x = [0 0 0 0 0 0]';
tau_w = [0 0 0]';
B = [1 1; 0.26 -0.26];
tau = [150;4];

for k = 1:Ns
    t = (k-1)*ts;
    
    T = B\tau;
    
    [x, Tf] = USV( x, T, tau_w, ts, t, tf );
    
    xout(k,:) = [t x' T' Tf' ];
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


figure(1);hold on
xrange=[-10 3 20]; yrange = [-10 3 20];
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
set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
axis([Xmin Xmax,Ymin Ymax]);
grid
hold off;

figure(2)
plot(t,u,'linewid',2);